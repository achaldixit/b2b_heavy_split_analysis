
# instead of doing this pass truck_speed during the CapClustring call
# Fix later
truck_speed <- function(cnid){
  #
  #
  #
  #
  return(truck_speed)
}

CapClusteringWithTSP <-
  function(lm_sc_adjusted_cluster,
           truck_speed = 9,
           dates,
           #the_date = min(lm_sc_adjusted_cluster$date, na.rm = T),
           #end_date = max(lm_sc_adjusted_cluster$date, na.rm = T),
           capacity_constraint = 1000) {
    library(dplyr, warn.conflicts = FALSE)
    time_per_drop <- 20 #mins
    
    time_per_kg <- 0.036 #mins
    
    # Suppress summarise info
    options(dplyr.summarise.inform = FALSE)
    
    #This is an approximate algo, so setting the capacity constraint as 1000kgs.
    capacity_constraint <- 1000
    
    full_trip_details <- data.frame()
    
    # the_date <- min()
    # end_date <- the_date + 1
    #end_date <- max(lm_selected_sc$date, na.rm=T)
    
    #Run the loop for all the dates.
    
    
    # Handle the specific date condition
    for(input_date in dates) {
      
      input_date <- as.Date(input_date, origin='1970-01-01')
      
      print(input_date)
      
      lm_sc_adjusted_cluster %>%
        filter(date == input_date) -> daily_fmlm
      
      #Get the bigger clusters for that day
      
      daily_clusters <- unique(daily_fmlm$loc_cluster)
      
      #For each of these bigger clusters, the capacitated clustering will be done
      #i_loc_cluster == iteration counter
      
      for (i_loc_cluster in unique(daily_fmlm$loc_cluster)) {
        picked_cluster_details <- data.frame()
        
        #The bigger cluster that is selected for the day.
        daily_fmlm %>%
          filter(loc_cluster == i_loc_cluster) -> daily_selected_cluster
        
        #These points will have a separate dedicated truck as their load is more than the capacity constraints
        
        daily_selected_cluster %>%
          filter(demand > capacity_constraint) -> dedicated_trucks
        
        #Rest of the shipments are ordered according to their demand. This arranging is a crucial part of the algo.
        
        daily_selected_cluster %>%
          filter(demand < capacity_constraint) %>%
          arrange(desc(demand)) %>%
          as.data.frame() %>%
          rowid_to_column(var = 'id') -> dropping_points_ordered
        
        #The capacitated clusters will only form when there are enough points to form the clusters with the subsequent loads.
        
        #If the total demand of the rest of the dropping points is less than 300 kgs, then those points will be added with the dedicated trucks.
        
        if (sum(dropping_points_ordered$demand) < 300) {
          dropping_points_ordered %>%
            bind_rows(dedicated_trucks) %>%
            select(-id) %>%
            rowid_to_column(var = 'id') %>%
            mutate(
              cluster = 1,
              cen_lat = mean(lat),
              cen_long = mean(long)
            ) -> capacitated_clustered_trips
          
          #Now if the demand is more than 300 but less than the capacity constraint then they can be clubbed together directly to form a cluster.
        } else if (sum(dropping_points_ordered$demand) > 300 &
                   sum(dropping_points_ordered$demand) < (capacity_constraint + 50)) {
          dropping_points_ordered %>%
            select(-id) %>%
            mutate(
              cluster = 1,
              cen_lat = mean(lat),
              cen_long = mean(long)
            ) %>%
            bind_rows(dedicated_trucks %>%
                        mutate(
                          cluster = 2,
                          cen_lat = lat,
                          cen_long = long
                        )) %>%
            rowid_to_column(var = 'id') -> capacitated_clustered_trips
        } else{
          if (nrow(dropping_points_ordered) > 0) {
            #Approximate number of clusters that will be formed. This is also the total number of trips.
            
            number_of_clusters <-
              ceiling(sum(dropping_points_ordered$demand) / (capacity_constraint * 0.8)) +
              2
            
            #For the first iteration, the centroids will be the points with larger loads
            
            dropping_points_ordered %>%
              slice(1:number_of_clusters) %>%
              select(lat, long) %>%
              rowid_to_column(var = 'cluster') -> centroids
            
            #The if condition here basically states that if the number of centroids and the number of points are equal, then there is no need to run the algo.
            #Basically, we can then treat each of the points as the center of the cluster.
            
            if (nrow(centroids) != nrow(dropping_points_ordered)) {
              n_iteration <- 1
              
              while (n_iteration <= 20) {
                #The cluster capacity = Max capacity of a cluster - The capacity already assigned
                
                cluster_capacity <-
                  capacity_constraint - as.vector(dropping_points_ordered[1:number_of_clusters, 'demand'])
                
                #All the points will have 'NA' from the beginning.
                
                cluster_list <-
                  rep(NA, nrow(dropping_points_ordered))
                
                #The centroids will have a defined clusters from the beginning. So the NAs for them will be replaced by numbers. These are the unique identification field for clusters.
                
                cluster_list[1:nrow(centroids)] <- 1:nrow(centroids)
                
                #Rest will be used for the capacitated clustering.
                
                dropping_points_ordered %>%
                  slice(c(number_of_clusters + 1:n())) -> points_for_clustering
                
                #i <- the iteration counter. Highlights the number of points for which we will run the algo.
                #The initial number is equal to the number of centroids. And the maximum value it can take is equal to the number of points.
                
                i <- nrow(centroids)
                
                while (i <= (nrow(dropping_points_ordered) - 1)) {
                  #j is the index that will be used here. It signifies the index for nearest cluster.
                  #It will start from 2. Suppose we are computing the distance matrix of point p1 and n-1 other points. So we will have a matrix of 1 * n including point p1.
                  #Now the lowest distance will always be 0 because the matrix will have p1 point in it. So the index of the 2nd nearest point will be 2.
                  
                  j <- 2
                  
                  #Randomly select one point.
                  
                  selected_point <-
                    slice_sample(points_for_clustering)
                  
                  #Select the unique id of that pont
                  
                  selected_point_id <- selected_point$id
                  
                  #Remove that point from the potential group of points being considered for clustering.
                  
                  points_for_clustering %>%
                    filter(id != selected_point_id) -> points_for_clustering
                  
                  #Now basically we will assign the randomly selected point to the nearest cluster considering the cluster has the capacity.
                  
                  selected_point_demand <- selected_point$demand
                  
                  selected_point_location <-
                    selected_point[, c('lat', 'long')]
                  
                  dist_from_clusters <-
                    geodist(x = as.matrix(selected_point_location),
                            y = as.matrix(centroids)) / 1000
                  
                  priority <-
                    dist_from_clusters / selected_point_demand
                  
                  nearest_cluster <- which.min(priority)
                  
                  nearest_cluster_capacity <-
                    cluster_capacity[nearest_cluster]
                  
                  if (selected_point_demand < nearest_cluster_capacity) {
                    cluster_list[selected_point_id] <- nearest_cluster
                    
                    dropping_points_ordered$cluster <- cluster_list
                    
                    cluster_capacity[nearest_cluster] <-
                      nearest_cluster_capacity - selected_point_demand
                    
                  } else{
                    #If the cluster doesn't have the capacity, then assign it to the 2nd nearest cluster.
                    #Run the loop until the point is assigned.
                    
                    while (is.na(cluster_list[selected_point_id])) {
                      nearest_cluster <- order(unlist(priority))[j]
                      nearest_cluster_capacity <-
                        cluster_capacity[nearest_cluster]
                      
                      if (selected_point_demand < nearest_cluster_capacity) {
                        cluster_list[selected_point_id] <- nearest_cluster
                        
                        dropping_points_ordered$cluster <-
                          cluster_list
                        
                        cluster_capacity[nearest_cluster] <-
                          nearest_cluster_capacity - selected_point_demand
                      } else{
                        j <- j + 1
                      }
                      
                      
                    }
                  }
                  
                  
                  i <- i + 1
                }
                
                
                dropping_points_ordered %>%
                  group_by(cluster) %>%
                  summarise(lat = mean(lat),
                            long = mean(long)) -> centroids
                
                n_iteration <- n_iteration + 1
              }
              
              
              picked_cluster_details %>%
                bind_rows(
                  centroids %>%
                    rename(cen_lat = lat, cen_long = long) %>%
                    left_join(dropping_points_ordered, by = 'cluster')
                ) -> picked_cluster_details
              
            } else{
              dropping_points_ordered %>%
                mutate(
                  cluster = id,
                  cen_lat = lat,
                  cen_long = long
                ) -> picked_cluster_details
            }
            
            capacitated_clustered_trips <- data.frame()
            
            for (x_date in unique(picked_cluster_details$date)) {
              picked_cluster_details %>%
                filter(date == x_date) -> temp_regular_cluster
              
              dedicated_trucks %>%
                filter(date == x_date) %>%
                rowid_to_column(var = 'id') %>%
                mutate(
                  id = max(temp_regular_cluster$cluster) + id,
                  cluster = id,
                  cen_lat = lat,
                  cen_long = long
                ) -> temp_truck
              
              
              capacitated_clustered_trips %>%
                bind_rows(temp_regular_cluster, temp_truck) -> capacitated_clustered_trips
            }
            
            #Further modification
            
            #Select the bad trips/trucks, where the load is less than 400
            
            capacitated_clustered_trips %>%
              group_by(cluster) %>%
              summarise(demand = sum(demand), .groups = 'drop') %>%
              filter(demand < 400) %>%
              pull(cluster) -> bad_trucks
            
            capacitated_clustered_trips %>%
              filter(cluster %in% bad_trucks) -> bad_truck_details
            
            #Select the good trips/trucks
            
            capacitated_clustered_trips %>%
              group_by(cluster) %>%
              summarise(demand = sum(demand), .groups = 'drop') %>%
              filter(demand > 400) %>%
              pull(cluster) -> good_trucks
            
            capacitated_clustered_trips %>%
              filter(cluster %in% good_trucks) -> good_truck_details
            
            bad_trucks <- c(bad_trucks)
            
            modified_bad_truck <- data.frame()
            
            #This code chunk will only be executed if there is at-least one bad truck.
            #The load of the bad cluster will be allocated to the nearest good cluster.
            
            if (nrow(bad_truck_details) > 0) {
              if (nrow(good_truck_details) == 0) {
                bad_truck_details %>%
                  mutate(cluster = 1) %>%
                  group_by(cluster) %>%
                  mutate(cen_lat = mean(lat),
                         cen_long = mean(long)) %>%
                  ungroup() -> capacitated_clustered_trips
                
                
              } else{
                for (i_bad_trucks in 1:nrow(bad_truck_details)) {
                  bad_truck_details[i_bad_trucks, ] -> selected_bad_truck
                  
                  truck_distances <-
                    rdist(x1 = selected_bad_truck[, c('lat', 'long')], x2 = good_truck_details[, c('cen_lat', 'cen_long')]) *
                    110
                  
                  nearest_truck <- which.min(truck_distances)
                  
                  good_truck_details %>%
                    slice(nearest_truck) -> nearest_good_truck
                  
                  selected_bad_truck$cluster <-
                    nearest_good_truck$cluster
                  
                  modified_bad_truck %>%
                    bind_rows(selected_bad_truck) -> modified_bad_truck
                }
                
                good_truck_details %>%
                  bind_rows(modified_bad_truck) %>%
                  group_by(cluster) %>%
                  mutate(cen_lat = mean(lat),
                         cen_long = mean(long)) %>%
                  ungroup() -> capacitated_clustered_trips
                
              }
              
              
            }
            
          } else{
            capacitated_clustered_trips <- data.frame()
            
            for (x_date in unique(dedicated_trucks$date)) {
              dedicated_trucks %>%
                filter(date == x_date) %>%
                rowid_to_column(var = 'id') %>%
                mutate(
                  cluster = id,
                  cen_lat = lat,
                  cen_long = long
                ) -> dedicated_trucks
              
              capacitated_clustered_trips %>%
                bind_rows(dedicated_trucks) -> capacitated_clustered_trips
            }
          }
          
        }
        
        capacitated_clustered_trips %>%
          group_by(cluster, km_cen_lat, km_cen_long) %>%
          count() %>%
          ungroup() %>%
          group_by(cluster) %>%
          mutate(perc = n / sum(n)) %>%
          slice_max(perc, n = 1) %>%
          ungroup() %>%
          distinct(cluster, km_cen_lat, km_cen_long) %>%
          left_join(capacitated_clustered_trips %>% select(-km_cen_lat,-km_cen_long),
                    by = "cluster") -> capacitated_clustered_trips
        
        #TSP for each trip to calculate the time to serve the cluster.
        
        picked_cluster_details <- capacitated_clustered_trips
        picked_tsp <- data.frame()
        
        for (x_cluster in unique(picked_cluster_details$cluster)) {
          picked_cluster_details %>%
            filter(cluster == x_cluster) -> selected_cluster
          
          if (nrow(selected_cluster) > 1) {
            distance_list <-
              osrmTable(loc = selected_cluster[, c('id', 'long', 'lat')], measure = 'distance')
            
            tsp_matrix <-
              ((distance_list$distances) / 1000) * (60 / truck_speed(cnid_input))
            
          } else{
            tsp_matrix <- as.matrix(0, 1, 1)
          }
          
          atsp <- ATSP(tsp_matrix)
          tour <- solve_TSP(atsp, method = "nn")
          
          picked_tsp %>%
            bind_rows(
              data.frame(
                date = as.Date(input_date),
                cluster = x_cluster,
                tsp_distance = (tour_length(tour) * (truck_speed(cnid_input) / 60)),
                tsp_time = tour_length(tour)
              )
            ) -> picked_tsp
          
        }
        
        
        picked_cluster_details %>%
          group_by(cluster,
                   date,
                   cen_lat,
                   cen_long,
                   purpose,
                   loc_cluster,
                   km_cen_lat,
                   km_cen_long) %>%
          summarise(
            count = n(),
            demand = sum(demand),
            .groups = 'drop'
          ) %>%
          ungroup() %>%
          inner_join(picked_tsp, by = c("cluster", "date")) %>%
          mutate(serving_time = ((count * time_per_drop) + (demand * time_per_kg))) -> picked_full_details
        
        full_trip_details %>%
          bind_rows(picked_full_details) -> full_trip_details
      }
      
      print(paste0('Done for the day: ', input_date, ' ', Sys.time()))
      
      #Move on to the next day.
      
    }
    
    return(full_trip_details)
  }
