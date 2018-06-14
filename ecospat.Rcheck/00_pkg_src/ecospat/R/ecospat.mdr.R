## Written by Olivier Broennimann. Departement of Ecology and Evolution (DEE).
## University of Lausanne. Switzerland. May 2011.  
## DESCRIPTION
## mdr(data,xcol,ycol,datecol,mode,rep,mean.date.error,fixed.sources.rows)
## mdr (minimum dispersal routes) is a function that implement a minimum cost
## arborescence approache to analyse the invasion routes of a species from dates
## occurrence data.  The function draw an incoming route to each occurence from the
## closest occurrence already occupied (with a previous date) The function allow to
## substract a random number of time (years) to the observed dates from a truncated
## negative exponential distribution. It is possible to run several iterations and to
## get boostrap support for each route. A list is returned by the function with in
## positon [[1]], a datafame with each row corresponding to a route (with new/old
## coordinates, new/old dates, length of the route, timespan, dispersal rate), in
## position [[2]] the total route length, in position [[3]] the median dispersal rate
## and in position [[4]] the number of outgoing nodes (index of clustering of the
## network) itexp and rtexp functions are small functions to set a truncated negative
## exponential distribution. The approach is decribed in HOrdijk & Broennimann 2012 J
## Theor Biol data = dataframe with occurence data. Each row correspond to an
## occurrence. xcol = the column in data containing x coordinates ycol = the column in
## data containing y coordinates datecol = the column in data containing dates mode =
## 'observed', 'min' or 'random'. 'observed' calculate routes using real dates. 'min'
## reorder dates so the the total length of the routes are minimal. 'random'
## reatribute dates randomly. rep = number of iteration of the analyse. if > 1,
## boostrap support for each route is provided mean.date.error = mean number of years
## to substract to observed dates. More precisely, it is the mean of the truncated
## negative exponential distribution from which the time to be substracted is randomly
## sampled.  fixed.sources.rows = the rows in data (as a vector) corresponding to
## source occurrence(s) that initiated the invasion(s). No incoming routes are not
## calculated for sources.



ecospat.mdr <- function(data, xcol, ycol, datecol, mode, rep, mean.date.error, fixed.sources.rows) {
  
  # symetrical matrix storing how many time a route is found between two occurrences
  links <- matrix(rep(0, nrow(data)^2), nrow = nrow(data), ncol = nrow(data), dimnames = list(c(1:nrow(data)), 
    c(1:nrow(data))))
  # symetrical matrix storing the sum of timespan of dispersal event for a route between
  # two occurrences
  timespan <- matrix(rep(0, nrow(data)^2), nrow = nrow(data), ncol = nrow(data), dimnames = list(c(1:nrow(data)), 
    c(1:nrow(data))))
  itexp <- function(u, m, t) {
    -log(1 - u * (1 - exp(-t * m)))/m
  }  #where u is the quantile, m is the rate, and t is the level of truncation
  rtexp <- function(n, m, t) {
    itexp(runif(n), m, t)
  }
  for (iter in 1:rep) {
    
    train <- data[, c(xcol, ycol, datecol)]  #data to be used
    if (mean.date.error == 0) 
      {
        e <- rep(0, nrow(train))
      }  # be aware that if too neighbor occurrences have the same date, a route is not possible.
    # you may have to set mean.date.error to a very small value
    if (mean.date.error > 0) {
      b <- rexp(nrow(train), 1/mean.date.error)
      q95 <- quantile(b, 0.95)
      
      
      e <- rtexp(nrow(train), 1/mean.date.error, q95)  # random sample of errors in a truncated (<0.95 quantile) neg. exp. distrib.
    }
    s <- fixed.sources.rows
    dist <- as.matrix(dist(train[, 1:2]))  # matrix a euclidian distances between sites
    
    if (mode == "observed") 
      {
        train[, 3] <- train[, 3] - e
      }  # substract the errors from observed dates
    if (mode == "random") 
      {
        train[, 3] <- sample(train[, 3] - e, nrow(train))
      }  # substract the errors from randomized dates
    if (mode == "min") {
      d <- dist
      d[d == 0] <- 99999
      if (is.null(s)) {
        ini <- sample(1:nrow(train), 1)
      }
      if (!is.null(s)) {
        ini <- s
      }
      while (length(ini) < nrow(train)) {
        min <- min(d[-ini, ini])
        a <- which(d[, ini] == min, arr.ind = T)[1]
        ini <- c(ini, a)
      }
      train[, 3] <- sort(train[, 3])[order(ini)] - e  # substract the errors from sorted dates
    }
    
    if (is.null(s)) {
      rows <- 1:nrow(train)
    }
    if (!is.null(s)) {
      rows <- (1:nrow(train))[-s]
    }
    
    for (i in rows) {
      prevsites <- which(train[, 3] < train[i, 3])  # find sites with a previous date
      
      if (length(prevsites) == 1) {
        links[i, prevsites] <- links[i, prevsites] + 1  # store the route
        timespan[i, prevsites] <- timespan[i, prevsites] + train[i, 3] - train[prevsites, 
          3]  # store the timespan of route
      }
      if (length(prevsites) > 1) {
        closestprevsite <- which(row.names(train) == as.integer(names(which.min(dist[i, 
          prevsites])))[1])  # find the closest site among those with previous dates
        links[i, closestprevsite] <- links[i, closestprevsite] + 1  # store the route
        timespan[i, closestprevsite] <- timespan[i, closestprevsite] + train[i, 
          3] - train[closestprevsite, 3]  # store the timespan of route
      }
    }
  }
  
  routes <- data.frame(matrix(nrow = sum(links > 0), ncol = 10, dimnames = list(c(), 
    c("Xold", "Yold", "Xnew", "Ynew", "dist", "dateold", "datenew", "timespan", "dispersal.rate", 
      "bootstrap.value"))))
  
  m = 1  #m=each individual route
  for (k in 1:nrow(train)) {
    for (l in 1:nrow(train)) {
      if (links[l, k] > 0) {
        routes[m, 1:2] <- train[k, 1:2]
        routes[m, 3:4] <- train[l, 1:2]
        routes[m, 5] <- sqrt((routes[m, 1] - routes[m, 3])^2 + (routes[m, 2] - 
          routes[m, 4])^2)
        routes[m, 6] <- train[k, 3]
        routes[m, 7] <- train[l, 3]
        routes[m, 8] <- timespan[l, k]/rep
        routes[m, 9] <- routes[m, 5]/routes[m, 8]
        routes[m, 10] <- links[l, k]/rep
        m = m + 1
      }
    }
  }
  
  tot.routes.length <- sum(routes$dist * routes$bootstrap.value)
  disp.rate <- median(rep(routes$dispersal.rate, routes$bootstrap.value * 100))
  out.nodes <- nrow(unique(routes[, 1:2]))
  
  routes <- routes[order(routes[, 10]), ]
  
  return(list(routes, tot.routes.length, disp.rate, out.nodes))
}