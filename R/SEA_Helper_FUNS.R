## Demo script for R burnr package
## Chris Guiterman
## chguiterman@email.arizona.edu
## 2018-5-27

## Some useful functions for SEA graphics
require(dplyr)

n_events <- function(df){
  n.val <- nrow(df)
  out <- substitute(italic(n)~"="~n.val,
                    list(n.val = format(n.val)))
  as.character(as.expression(out))
}

sea_sig <- function(x){
  out <- x %>% mutate(sig = if_else(mean <= lower_95_perc, "sig_neg",
                                    if_else(mean > lower_95_perc & mean <=0, "neg",
                                            if_else(mean > 0 & mean < upper_95_perc, "pos", "sig_pos"))))
  out$sig <- factor(out$sig, levels=c("sig_neg", "neg", "pos", "sig_pos"))
  return(out)
}

sea_cols <- c('#ca0020', '#f4a582', '#92c5de', '#0571b0')
names(sea_cols) <- c("sig_neg", "neg", "pos", "sig_pos")