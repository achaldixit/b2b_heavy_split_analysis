library(tidyverse)
library(data.table) 
library(lubridate)
library(geodist)
library(ompr)
library(ompr.roi)
library(ROI.plugin.symphony)
library(ggrepel)
# Mention parameters in chunks to get notebook outputs
library(palmerpenguins)
library(TSP)
library(fields)
library(Rglpk)
library(ROI.plugin.glpk)
library(osrm)

# make sure to setwd() and keep the file in same directory
source("/home/achal/R_codes/b2b_heavy_split_analysis/capacitated_clustering_tsp.R")

PATH <-
  "/rstudio-dwh-s3/shibaprasad/sc_pickup_drop_august2022.RData"

PATH_OUTPUT = "/home/shibaprasad/b2b_heavy_split/analysis_output/"

data <- get(load(PATH))

#cnid_input <- 'INTNAQRV'

fmlm_sc %>%
  inner_join(fmlm_sc %>%
               distinct(cnid, product_type) %>%
               count(cnid) %>%
               filter(n == 2),
             by = 'cnid') %>%
  filter(!str_detect(ptype, 'FC')) %>%
  group_by(cn, cnid) %>%
  summarise(cn_demand = sum(weight, na.rm = T), .groups = 'drop') %>%
  arrange(desc(cn_demand)) %>%
  mutate(cum_demand_perc = 100 * cumsum(cn_demand) / sum(cn_demand)) %>%
  filter(cum_demand_perc < 81) %>%
  pull(cnid) ->
  cnid_list

#dput(cnid_list[cnid_list!=unique(analysis_output$cnid)])


c("INDLAAFL", "INTNAQRV", "INMHASHG", "INTNARDZ", "INDLAOSF","INWBAPMU", "INKAASBY",
  "INTSARDI", "IND700088AAA", "IND160102AAC","INDLAOUR", "INUPAREB", "INASAQXA", "INTNAREC", "INKAAASD", "INMHABET",
  "INBRAANC", "INHRAPQY", "INKAAQXS", "INKLARDT", "INHRATBM", "INUPARDS",
  "INWBAPEX", "INKAAQYO", "INHRAQTD", "INJHARDW", "IND560083AAB",
  "INDLAQPF", "INHRASJG", "INHRAAKG", "INTSAMWA", "INMHAQVB", "INMHAMVN",
  "INMHARNO", "INUKARXY", "INMHAQRF", "INMHAAXB", "INRJAOFY", "INKAANKD",
  "INDLARQO", "INAPAOWP", "INGJAQBZ", "INMHAAXV", "IND600044AAE",
  "INKAAQXV", "INMHATIR", "INKAAOHM", "INWBARGV", "INTSAOWM", "INMHAOIJ", 
  "INMPARLR", "INTSAMNQ", "INMHAQSZ", "INMHAAXI", "INGJAAOH", "INORAMZN", 
  "INUPAREG", "INPBAREH", "INMPARDO", "INTNAQRY", "INTNASTA", "INUPAOIC", 
  "INKAAQOA", "INMHAMVH", "INTSAQQK", "INMHAOMN", "INDLANOY", "INMHABFA", 
  "INHPAAEH", "INMHAQQO", "INTNAQUD", "INORARDF", "INUPAQQA", "INUPAAXU", 
  "IN700069A1A", "INTSASAI", "IND570020AAC", "INRJAREN", "INJKANMI", 
  "INORAAKD", "INMHAQZI", "INRJAROK", "INGJAAOD", "INGJASIZ", "INKAAQXU", 
  "INUPAPEI")-> cnid_list

# Scope of this function is to obtain trip cost for a pdt, cnid, for a given month.
# Assumes all of the data at source is utilized.

GetTripCosts <- function(data,
                         prod_type,
                         cnid_input,
                         dates,
                         capacity_constraint = 1000,
                         lowest_cluster_demand = 1000,
                         lowest_large_cluster_demand = 3000) {
  time_per_drop <- 20 #mins
  
  time_per_kg <- 0.036 #mins
  
  # load the data from path source
  fmlm_sc <- data
  
  fmlm_sc %>%
    mutate(date = as.Date(dispatch_date)) %>%
    filter(date %in% dates) %>%
    filter(cnid %in% cnid_input & purpose == 'LM') %>%
    filter(product_type %in% prod_type) %>%
    filter(!is.na(customer_lat) &
             !is.na(customer_long)) -> lm_selected_sc
  
  print(cnid_input)
  
  q1_lat <-
    median(lm_selected_sc$customer_lat) - 2 * sd(lm_selected_sc$customer_lat)
  q2_lat <-
    median(lm_selected_sc$customer_lat) + 2 * sd(lm_selected_sc$customer_lat)
  
  q1_long <-
    median(lm_selected_sc$customer_long) - 2 * sd(lm_selected_sc$customer_long)
  q2_long <-
    median(lm_selected_sc$customer_long) + 2 * sd(lm_selected_sc$customer_long)
  
  lm_selected_sc %>%
    filter((customer_lat > q1_lat) &
             (customer_lat < q2_lat) &
             (customer_long > q1_long) & (customer_long < q2_long)
    ) %>%
    group_by(cnid, date) %>%
    mutate(weight = ifelse(is.na(weight), median(weight, na.rm = TRUE), weight)) %>%
    ungroup() %>%
    group_by(customer_lat,
             customer_long,
             date,
             product_type,
             purpose,
             dispatch_id) %>%
    summarise(demand = sum(weight), .groups = 'drop') %>%
    rename(lat = customer_lat, long = customer_long) -> lm_selected_sc
  
  
  lm_selected_sc %>%
    group_by(date) %>%
    summarise(demand = sum(demand, na.rm = T), .groups = 'drop') %>%
    arrange((demand)) %>%
    pull(date) -> list_of_dates
  
  # Dispatch_ID gives shitty clustering with insane overlap
  # lm_selected_sc %>%
  #   mutate(loc_cluster=as.numeric(as.factor(dispatch_id))) %>%
  #   group_by(loc_cluster) %>%
  #   mutate(km_cen_lat=mean(lat),
  #          km_cen_long=mean(long)) %>%
  #   ungroup() %>%
  #   arrange(loc_cluster)-> lm_selected_sc
  
  # Kmeans CLustering
  set.seed(4)
  
  #Setting a lowest cluster demand for large cluster which will give us the number of clusters.
  
  d <- lm_selected_sc[,c('lat','long')]
  n_centers <- ceiling(sum(lm_selected_sc$demand, na.rm=T)/(30*lowest_large_cluster_demand))
  
  print(paste0("N is :",n_centers))
  
  # fmlm_sc %>%
  #   filter(cnid %in% cnid_input) %>%
  #   pull(dispatch_id) %>%
  #   n_distinct() -> n_centers
  
  kmeans <- kmeans(x=d, centers=round(n_centers))
  
  lm_selected_sc$loc_cluster <- kmeans$cluster
  
  lm_selected_sc %>%
    group_by(loc_cluster) %>%
    mutate(km_cen_lat=mean(lat),
           km_cen_long=mean(long)) %>%
    ungroup()-> lm_selected_sc
  
  lm_selected_sc
  
  lm_sc_adjusted_cluster <- data.frame()
  
  # Cluster merging
  for (the_date in unique(lm_selected_sc$date)) {
    the_date <- as.Date(the_date, origin = '1970-01-01')
    
    print(the_date)
    
    lm_selected_sc %>%
      filter(date == the_date) -> selected_date
    
    selected_date %>%
      group_by(loc_cluster, km_cen_lat, km_cen_long) %>%
      summarise(demand = sum(demand), .groups = 'drop') -> demand_cluster_date
    
    demand_cluster_date %>%
      filter(demand <= lowest_cluster_demand) -> poor_clusters
    
    if (nrow(poor_clusters) == 0) {
      selected_date %>%
        bind_rows(lm_sc_adjusted_cluster) -> lm_sc_adjusted_cluster
      next
    }
    
    demand_cluster_date %>%
      filter(demand > lowest_cluster_demand) -> better_clusters
    
    if (nrow(better_clusters) == 0) {
      selected_date %>%
        mutate(loc_cluster = 1) %>%
        bind_rows(lm_sc_adjusted_cluster) -> lm_sc_adjusted_cluster
      
    } else{
      for (i in 1:nrow(poor_clusters)) {
        nearest_cluster <-
          which.min(rdist(x1 = poor_clusters[i, c('km_cen_lat', 'km_cen_long')], x2 =
                            better_clusters[, c('km_cen_lat', 'km_cen_long')]) * 110)
        
        selected_poor_cluster <- poor_clusters$loc_cluster[i]
        
        better_clusters %>%
          slice(nearest_cluster) -> assigned_better_cluster
        
        selected_date$km_cen_lat[selected_date['loc_cluster'] == selected_poor_cluster] <-
          assigned_better_cluster$km_cen_lat
        selected_date$km_cen_long[selected_date['loc_cluster'] == selected_poor_cluster] <-
          assigned_better_cluster$km_cen_long
        
        selected_date %>%
          mutate(
            loc_cluster = ifelse(
              loc_cluster == selected_poor_cluster,
              assigned_better_cluster$loc_cluster,
              loc_cluster
            )
          ) -> selected_date
        
      }
      
      selected_date %>%
        bind_rows(lm_sc_adjusted_cluster) -> lm_sc_adjusted_cluster
      
    }
  }
  
  # Fix later: Include the reference table for truck speed
  GetTruckSpeed <- function(cnid_input){

    truck_speeds = list("IN700069A1A" = 5)
    
    return(truck_speeds[cnid_input][[1]])
  }
  
  truck_speed <- GetTruckSpeed(cnid_input)
  print("Truck Speed")
  print(truck_speed)
  # CapacitatedClusteringWithTSP source the R file for this
  
  # lm_sc_adjusted_cluster %>%
  #   filter(date >= as.Date('2022-08-11') &
  #            date <= as.Date('2022-08-12')) -> d
  
  CapClusteringWithTSP(lm_sc_adjusted_cluster, 
                       truck_speed,
                       dates = dates
                       ) -> full_trip_details
  
  # add SC Location
  
  full_trip_details %>%
    mutate(
      plat = fmlm_sc %>% filter(cnid == cnid_input) %>% distinct(plat) %>% pull(),
      plong = fmlm_sc %>% filter(cnid == cnid_input) %>% distinct(plong) %>% pull(),
      stem_time = 2 * (abs(plat - cen_lat) + abs(plong - cen_long)) *
        110 * 60 / truck_speed,
      cluster_time = stem_time + serving_time + tsp_time
    ) -> full_trip_details
  
  # Allocate Trucks as per the load
  full_trip_details %>%
    mutate(
      load_category = case_when(
        demand <= 1000 ~ 'ace',
        demand > 1000 & demand <= 2500 ~ 'truck12ft',
        demand > 2500 & demand <= 3500 ~ 'truck14ft',
        demand > 3500 & demand <= 5500 ~ 'truck17ft',
        demand > 5500 & demand <= 8500 ~ 'largetruck20ft',
        demand > 8500 & demand <= 18000 ~ 'largetruck24ft',
        demand > 18000 & demand <= 24500 ~ 'largetruck40ft',
        demand > 24500 ~ 'truck46ft'
      )
    ) -> full_trip_details
  
  vehicle_cost <-
    data.frame(
      vehicle = c(
        'ace',
        'truck12ft' ,
        'truck14ft',
        'truck17ft',
        'largetruck20ft',
        'largetruck24ft',
        'largetruck40ft',
        'truck46ft'
      ),
      cost = c(2800, 3300, 3600, 4920, 6200, 10616, 13700, 14000)
    )
  
  # Knapsack the truck loads
  
  full_trip_cost <- data.frame()
  
  full_trip_details %>%
    left_join(vehicle_cost %>% rename(load_category = vehicle), by = "load_category") %>%
    mutate(
      contracted_shift = cluster_time / 720,
      trip_cost = ifelse(contracted_shift >= 1, contracted_shift * cost, cost)
    ) %>%
    mutate(cnid = cnid_input)-> full_trip_cost
  
  return(full_trip_cost)
}

# Function to run bulk analysis for Multiple SC's
RunProductAnalysis <- function(data, cnid_df, product_type) {
  analysis_output <- tibble()
  cnid_list <- unique(cnid_df$cnid)
  print(cnid_list)
  for (cnid_in in cnid_list) {
    
    # Obtain the dates to run for
    cnid_df %>%
      filter(cnid == cnid_in) ->
      cnid_single
      
    # Only run for the dates specified
    dates <- cnid_single$date
    cnid_output <-
      GetTripCosts(data, 
                   cnid_input = cnid_in, 
                   prod_type = product_type,
                   dates = dates) %>%
      mutate(product_type = toString(product_type))
    
    analysis_output %>%
      bind_rows(cnid_output) ->
      analysis_output
    
    Sys.umask(002)
    fname = str_glue(PATH_OUTPUT,"split_analysis_999",toString(product_type),"_", ".RData")
    save(analysis_output, file=fname)
    
    gc()
    
  }
  return(analysis_output)
}

load("/rstudio-dwh-s3/shibaprasad/b2b_heavy_split/trip_cost_comparison.RData")

trip_cost_comparison %>%
  mutate(dow = weekdays(date)) %>%
  filter(dow != 'Sunday' &
           date != as.Date('2022-08-15') & cost_diff_perc >= 0) %>%
  distinct(cn, cnid, date) ->
  sample_cnid

sample_cnid %>%
  filter(cnid == "IN700069A1A") ->
  sample

product_type = c("B2B")
#cnid_list = c('IN700069A1A')#, "INTNAQRV", "INMHASHG") #"INTNARDZ", "INDLAOSF")
#analysis_output <- RunProductAnalysis(data, cnid_list, product_type = product_type)

analysis_output <- RunProductAnalysis(data, sample, product_type = product_type)

Sys.umask(200)
product_type = str_replace(toString(product_type), ",","")
product_type = str_replace(product_type, " ","_")
fname = str_glue(PATH_OUTPUT,"split_analysis_",toString(product_type),"_",str_replace(toString(Sys.time())," ", "_"), ".RData")
save(analysis_output, file = fname)