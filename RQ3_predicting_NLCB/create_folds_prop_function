create_folds_prop <- function (data, V = 5, id_var = "id", seed = 123L, t0 = 4) 
{
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1L)
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
  set.seed(seed)
  
  # before t0
  df1 <- as.data.frame(data)
  df1 <- df1[df1$time_obs <= t0,]
  ids <- df1[[id_var]]
  unq_ids <- unique(ids)
  n <- length(unq_ids)
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  training1 <- testing1 <- vector("list", V)
  for (i in seq_along(training1)) {
    ind <- ids %in% unq_ids[splits[[i]]]
    training1[[i]] <- df1[!ind, ]
    testing1[[i]] <- df1[ind, ]
  }
  
  # after t0
  df2 <- as.data.frame(data)
  df2 <- df2[df2$time_obs > t0,]
  ids <- df2[[id_var]]
  unq_ids <- unique(ids)
  n <- length(unq_ids)
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  training2 <- testing2 <- vector("list", V)
  for (i in seq_along(training2)) {
    ind <- ids %in% unq_ids[splits[[i]]]
    training2[[i]] <- df2[!ind, ]
    testing2[[i]] <- df2[ind, ]
  }
  
  # combine datasets
  training <- testing <- vector("list", V)
  for (i in seq_along(training)) {
    training[[i]] <- rbind(training1[[i]], training2[[i]])
    testing[[i]] <- rbind(testing1[[i]], testing2[[i]])
  }
  
  # output
  list(training = training, testing = testing)
}
