B_Corrector <- function(i_vector)
{
    n_noise_sd <- 1.5   #cutoff threshold expressed as times of noise standard deviation, default = 1.1
    sd_interval <- 1    #window size is 41 (2*20+1)
    involve_interval <- 1   #3
    involve_interval_2 <- 3 #7
    average_interval <- 10  #60
    tmp_list <- list()
    #---------------------------------------------------------------------------
    #determination of noise standard deviation
    tmp_sd <- rep(0, length(i_vector))
    for (j in (sd_interval + 1):(length(i_vector) - sd_interval))
    {
        tmp_sd[j] <- sd(i_vector[(j - sd_interval):(j + sd_interval)])
    }
    for (j in 1:sd_interval)
    {
        tmp_sd[j] <- sd(i_vector[1:(j + sd_interval)])
    }
    for (j in (length(i_vector) - sd_interval + 1):length(i_vector))
    {
        tmp_sd[j] <- sd(i_vector[(j - sd_interval):length(i_vector)])
    }
    tmp_sd_no_zero <- tmp_sd[!tmp_sd == 0]
    tmp_median <- median(tmp_sd_no_zero)
    tmp_median1 <- tmp_median
    repeat
    {
        tmp_median2 <- median(tmp_sd_no_zero[tmp_sd_no_zero < (tmp_median1*2)])
        if (abs(tmp_median2 - tmp_median1 ) < (tmp_median1/1000))
        {
            break
        }   else
        {
            tmp_median1 <- tmp_median2
        }
    }
    noise_sd <- tmp_median2    #noise_sd: estimated noise standard deviation

    #---------------------------------------------------------------------------
    #classify signal/noise
    tmp_sd_mark <- rep(0, length(i_vector))

    for (j in (involve_interval + 1):(length(i_vector) - involve_interval))
    {
        if (tmp_sd[j] > (n_noise_sd*noise_sd))
        {
            tmp_sd_mark[(j - involve_interval):(j + involve_interval)] <- 1
        }
    }
    for (j in 1:involve_interval)
    {
        if (tmp_sd[j] > (n_noise_sd*noise_sd))
        {
            tmp_sd_mark[1:(j + involve_interval)] <- 1
        }
    }
    for (j in (length(i_vector) - involve_interval + 1):length(i_vector))
    {
        if (tmp_sd[j] > (n_noise_sd*noise_sd))
        {
            tmp_sd_mark[(j - involve_interval):length(i_vector)] <- 1
        }
    }
    #-------------------------------------------------------------------------------
    #Replace beginning and ending signal with nearby noise value.
    i_vector_h_t <- i_vector
    if (sum(tmp_sd_mark[1:(involve_interval_2 + 1)])<(involve_interval_2 + 1))
    {
          k <- 0
          tmp_sum <- 0
          for (i in 1:(involve_interval_2 + 1))
          {
               if (tmp_sd_mark[i] == 0)
               {
                  k <- k + 1
                  tmp_sum <- tmp_sum + i_vector[i]
               }
          }
          i_vector_h_t[1:(involve_interval_2 + 1)] <- tmp_sum/k
    }
    if (sum(tmp_sd_mark[(length(i_vector) - involve_interval_2):length(i_vector)]) < (involve_interval_2 + 1))
    {
          k <- 0
          tmp_sum <- 0
          for (i in (length(i_vector) - involve_interval_2):length(i_vector))
          {
               if (tmp_sd_mark[i] == 0)
               {
                  k <- k + 1
                  tmp_sum <- tmp_sum + i_vector[i]
               }
          }
          i_vector_h_t[(length(i_vector) - involve_interval_2):length(i_vector)] <- tmp_sum/k
    }


    if (sum(tmp_sd_mark[1:(involve_interval_2 + 1)])==(involve_interval_2 + 1))
    {
          i <- 1
          k <- 1
          tmp_sum <- 0
          while (k <= (2*involve_interval_2 + 1))
          {
               if (tmp_sd_mark[i] == 0)
               {
                  k <- k + 1
                  tmp_sum <- tmp_sum + i_vector[i]
               }
               i <- i + 1
          }
          i_vector_h_t[1:(involve_interval_2 + 1)] <- tmp_sum/(2*involve_interval_2 + 1)
    }
    if (sum(tmp_sd_mark[(length(i_vector) - involve_interval_2):length(i_vector)])==(involve_interval_2 + 1))
    {
          i <- length(i_vector)
          k <- 1
          tmp_sum <- 0
          while (k <= (2*involve_interval_2 + 1))
          {
               if (tmp_sd_mark[i] == 0)
               {
                  k <- k + 1
                  tmp_sum <- tmp_sum + i_vector[i]
               }
               i <- i - 1
          }
          i_vector_h_t[(length(i_vector) - involve_interval_2):length(i_vector)] <- tmp_sum/(2*involve_interval_2 + 1)
    }
    tmp_sd_mark[1:(involve_interval_2 + 1)] <- 0
    tmp_sd_mark[(length(i_vector) - involve_interval_2):length(i_vector)] <- 0
    #---------------------------------------------------------------------------
    #Replace signal segments with linear interpolation
    ptable_normal_noise <- i_vector_h_t
    counter <- 0
    for (j in (2+involve_interval_2):(length(i_vector)-involve_interval_2))
    {
        if (tmp_sd_mark[j]==1)
        {
          counter <- counter+1
        }
        if (tmp_sd_mark[(j-1)]==1 & tmp_sd_mark[j]==0)
        {
            if (counter==0)
            {
                next
            }   else
            {
                a <- mean(i_vector_h_t[(j-counter-1-involve_interval_2):(j-counter-1+involve_interval_2)])
                b <- mean(i_vector_h_t[(j-involve_interval_2):(j+involve_interval_2)])
                ptable_normal_noise[(j-counter):(j-1)] <- seq(a,b,length.out=counter)
                counter <- 0
            }
        }
    }

    #-------------------------------------------------------------------------------
    #Mean filter
    tmp_baseline <- rep(0, length(i_vector))
    for (i in 1:average_interval)
    {
        tmp_baseline[i] <- mean(ptable_normal_noise[1:(i + average_interval)])
    }
    for (i in (average_interval + 1):(length(i_vector) - average_interval))
    {
        tmp_baseline[i] <- mean(ptable_normal_noise[(i - average_interval):(i + average_interval)])
    }
    for (i in (length(i_vector) - average_interval + 1):length(i_vector))
    {
        tmp_baseline[i] <- mean(ptable_normal_noise[(i - average_interval):length(i_vector)])
    }
    #-------------------------------------------------------------------------------
    tmp_list <- list(bl = tmp_baseline, sp = i_vector - tmp_baseline, nsd = noise_sd)
    return(tmp_list)
}

