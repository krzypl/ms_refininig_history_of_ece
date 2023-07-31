#function for calculating relative contribution of diatioms grouped according to their salinity preference

diatom_eco_groups <- function(x, y) {
  dm <- rowSums(x[,cumsum(str_starts(colnames(x), "\\(")) == 1])
  db <- rowSums(x[,cumsum(str_starts(colnames(x), "\\(")) == 2])
  de <- rowSums(x[,cumsum(str_starts(colnames(x), "\\(")) == 3])
  df <- rowSums(x[,cumsum(str_starts(colnames(x), "\\(")) == 4])
  di <- rowSums(x[,cumsum(str_starts(colnames(x), "\\(")) == 5])
  z <- tibble(depth = y, M = dm, B = db, E = de, F = df, U = di)
}
