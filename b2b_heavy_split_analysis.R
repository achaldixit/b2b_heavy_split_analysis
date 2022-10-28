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
source("capacitated_clustering_tsp.R")

PATH_OUTPUT = "/home/achal/R_codes/b2b_heavy_split_analysis/analysis_output/"

PATH <-
  "/rstudio-dwh-s3/shibaprasad/sc_pickup_drop_august2022.RData"

data <- get(load(PATH))

cnid_input <- 'INTNAQRV'

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

# Scope of this function is to obtain trip cost for a pdt, cnid, for a given month.
# Assumes all of the data at source is utilized.

GetTripCosts <- function(data,
                         prod_type,
                         cnid_input,
                         capacity_constraint = 1000,
                         lowest_cluster_demand = 1000,
                         lowest_large_cluster_demand = 3000,
                         truck_speed = 9) {
  time_per_drop <- 20 #mins
  
  time_per_kg <- 0.036 #mins
  
  # load the data from path source
  fmlm_sc <- data
  
  fmlm_sc %>%
    mutate(date = as.Date(dispatch_date)) %>%
    filter(cnid %in% cnid_input & purpose == 'LM') %>%
    filter(product_type %in% prod_type) %>%
    filter(!is.na(customer_lat) &
             !is.na(customer_long)) -> lm_selected_sc
  
  
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
  
  n_centers <- floor(sum(lm_selected_sc$demand, na.rm=T)/(30*lowest_large_cluster_demand))
  
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
  
  
  # CapacitatedClusteringWithTSP source the R file for this
  
  # lm_sc_adjusted_cluster %>%
  #   filter(date >= as.Date('2022-08-11') &
  #            date <= as.Date('2022-08-12')) -> d
  
  CapClusteringWithTSP(lm_sc_adjusted_cluster, truck_speed) -> full_trip_details

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
        demand < 1000 ~ 'ace',
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
  
  for (i_date in unique(full_trip_details$date)) {
    multi_trip_cost <- data.frame()
    
    i_date <- as.Date(i_date, origin = '1970-01-01')
    
    full_trip_details %>%
      filter(date == i_date) -> fleets_per_day
    
    fleets_per_day %>%
      filter(cluster_time < 720) -> under_utilised_trips
    
    for (i_load_category in unique(under_utilised_trips$load_category)) {
      under_utilised_trips %>%
        filter(load_category == i_load_category) %>%
        rowid_to_column(var = 'i') -> under_utilised_fleets
      
      if (nrow(under_utilised_fleets) == 1) {
        under_utilised_fleets %>%
          inner_join(vehicle_cost %>%
                       rename(load_category = vehicle)) %>%
          mutate(vehicle_id = paste0(load_category, '_', date, '_', cnid_input, '_', 1)) %>%
          group_by(vehicle_id, cost) %>%
          summarise(total_time = sum(cluster_time),
                    .groups = 'drop') %>%
          mutate(
            contracted_shifts = total_time / 720,
            trip_cost = ifelse(contracted_shifts > 1, cost + (cost /
                                                                12), cost),
            cnid = cnid_input,
            date = i_date
          ) %>%
          bind_rows(multi_trip_cost) -> multi_trip_cost
        
      } else{
        fleet_total_time <- as.vector(under_utilised_fleets$cluster_time)
        
        vehicle_cost %>%
          filter(vehicle == i_load_category) %>%
          pull(cost) -> vehicle_rent_cost
        
        m <- nrow(under_utilised_fleets)
        n <- nrow(under_utilised_fleets)
        
        rent_model <- MIPModel() %>%
          add_variable(z[i, j],
                       i = 1:n,
                       j = 1:m,
                       type = "binary") %>%
          add_variable(v[j], j = 1:m, type = "binary") %>%
          set_objective(sum_expr(vehicle_rent_cost * v[j], j = 1:m), "min") %>%
          add_constraint(sum_expr(z[i, j], j = 1:m) == 1, i = 1:n) %>%
          add_constraint(z[i, j] <= v[j], i = 1:n, j = 1:m) %>%
          add_constraint(sum_expr(fleet_total_time[i] * z[i, j], i = 1:n) <= 780 * v[j], j = 1:m)
        
        
        rent_result <-
          solve_model(rent_model,
                      with_ROI("glpk", verbose = TRUE, tm_limit = 1000 * 10))
        
        # rent_result <-solve_model(rent_model,
        #           with_ROI(
        #             solver = "symphony",
        #             verbosity = -1,
        #             gap_limit = 10
        #           ))
        
        rent_result %>%
          get_solution(z[i, j]) %>%
          filter(value == 1) %>%
          select(i, j) %>%
          inner_join(under_utilised_fleets, by = 'i') %>%
          inner_join(vehicle_cost %>%
                       rename(load_category = vehicle), by = "load_category") %>%
          mutate(vehicle_id = paste0(load_category, '_', date, '_', cnid_input, '_', j)) %>%
          group_by(vehicle_id, cost) %>%
          summarise(total_time = sum(cluster_time),
                    .groups = 'drop') %>%
          mutate(
            contracted_shifts = total_time / 720,
            trip_cost = ifelse(contracted_shifts > 1, cost + (cost /
                                                                12), cost),
            cnid = cnid_input,
            date = i_date
          ) %>%
          bind_rows(multi_trip_cost) -> multi_trip_cost
      }
    }
    
    fleets_per_day %>%
      filter(cluster_time >= 720) %>%
      inner_join(vehicle_cost %>%
                   rename(load_category = vehicle), by = "load_category") %>%
      rowid_to_column(var = 'id') %>%
      mutate(id = id + nrow(multi_trip_cost)) %>%
      mutate(vehicle_id = paste0(load_category, '_', date, '_', cnid_input, '_', id)) %>%
      group_by(vehicle_id, cost) %>%
      summarize(total_time = sum(cluster_time), .groups = 'drop') %>%
      mutate(
        contracted_shifts = total_time / 720,
        trip_cost = (cost * contracted_shifts),
        date = i_date,
        cnid = cnid_input
      ) %>%
      bind_rows(multi_trip_cost) %>%
      bind_rows(full_trip_cost) -> full_trip_cost
    
    print(paste0('m-VRP simulation done for : ', i_date))
    
  }
  return(full_trip_cost)
}

# Function to run bulk analysis for Multiple SC's
RunProductAnalysis <- function(data, cnid_list, product_type) {
  analysis_output <- tibble()
  for (cnid in cnid_list) {
    cnid_output <-
      GetTripCosts(data, cnid_input = cnid, product_type = product_type) %>%
      mutate(product_type = toString(product_type))
    
    
    analysis_output %>%
      bind_rows(cnid_output) ->
      analysis_output
  }
  return(analysis_output)
}

product_type = c("B2B","Heavy")
cnid_list = c("INDLAAFL", "INTNAQRV", "INMHASHG") #"INTNARDZ", "INDLAOSF")
analysis_output <- RunProductAnalysis(data, cnid_list, product_type = product_type)

fname = str_glue(PATH_OUTPUT,"split_analysis_",toString(product_type),"_",str_replace(toString(Sys.time())," ", "_"), ".RData")
save(analysis_output, file = fname)