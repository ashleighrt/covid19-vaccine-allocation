#calculate bezier curve to describe changing infection risk over time 

bezier_fun<- function(params){
  
  y0 <- params$rr_0
  y1 <- params$rr_1
  y2 <- params$rr_2
  y3 <- params$rr_3
  
  p0 =  y0
  p1 = ( -5*y0 + 18*y1 -  9*y2 + 2*y3) / 6
  p2 = (  2*y0 -  9*y1 + 18*y2 - 5*y3) / 6
  p3 = y3
  p<-c(p0,p1,p2,p3)
  t<-seq(0,1, length.out=params$weeks_tot ) 
  
  y <- bezier::bezier(t,p) 
  y1 <- c(p0, y) #set incidence at week 0 = incidence at week 1
  return(y1)
}

vac_fun <- function(params){
  
  weeks_tot <- params$weeks_tot
  pop_size <- params$pop_size
  pct_new <- params$pct_new
  rel_protect_1_dose <- params$rel_protect_1_dose
  dose_gap <- params$dose_gap # delay between doses (weeks)
  steps <- params$steps
  start_weekly_doses <- params$start_weekly_doses
  end_weekly_doses <- params$end_weekly_doses
  start_pct_saved_fixed <- params$start_pct_saved_fixed
  start_pct_saved_switch <- params$start_pct_saved_switch
  waning_scalar <- params$waning_scalar
  delay <- params$delay
  ve <- params$ve
  scheme <- params$scheme
  base_inc <- params$base_inc
  switch_week <- params$switch_week
  
  # make functions for vaccine availability and allocation strategy
  sf_vac <- stepfun(x=params$switch_week, y = seq(params$start_weekly_doses, params$end_weekly_doses, 
                                                  length.out=2), right=FALSE)
  sf_allocate <- if(scheme=="fixed") rep(start_pct_saved_fixed, weeks_tot) else c(rep(start_pct_saved_switch, switch_week-1), rep(1-start_pct_saved_switch, switch_week-1), rep(start_pct_saved_fixed, weeks_tot-(switch_week-1)*2))
  rr_curve <- pull(rr_curves[params$rr_curve_num])

  num_vac <- rep(0, pop_size) # number of vaccine doses received
  week_vac <- rep(99999, pop_size) # week of vaccine receipt
  
  # store what state everyone is in 
  vac_pop_mat <- matrix(0, ncol = pop_size, nrow = weeks_tot +1 ) 
  week_first_vec <- rep(week_vac, pop_size) #week of first dose
  week_second_vec <- rep(week_vac, pop_size) #week of second dose
  
  # vectors to store number of vaccine doses available and saved
  new_vac_vec <- rep(0, weeks_tot+1)  # weekly new doses available
  new_reserve_vec <- rep(0, weeks_tot+1) # new vaccine saved for 2nd dose
  new_distributed_vec <- rep(0, weeks_tot+1) # vaccine given to previously unvaccinated people
  total_reserve_vec <- rep(0, weeks_tot+1) # total saved vaccines
  
  # initial states
  vac_pop_mat[1,] <- num_vac
  
  # for each time unit
  for(kk in 1:weeks_tot){
    # get number of doses available 
    new_vac <- sf_vac(kk) * pop_size # weekly total new doses available
    new_reserve_vac <- new_vac * (sf_allocate[kk]) # number of new vaccines set aside for second dose
    new_distributed_vac <- new_vac * (1 - sf_allocate[kk]) # number of new vaccines 
    total_reserve_vac <- total_reserve_vec[kk] + new_reserve_vac
    
    #allocate first doses to unvaccinated people
    vec <- vac_pop_mat[kk,] 
    if(sum(vec)==0){
      vec[1:new_distributed_vac] <- 1
      week_first_vec[1:new_distributed_vac] <- kk
    } else{
      week_first_vec[(max(which(vec>0))+1):(max(which(vec>0))+new_distributed_vac)] <- kk
      vec[(max(which(vec>0))+1):(max(which(vec>0))+new_distributed_vac)] <- 1
    }
    
    #allocate second doses to those with first dose
    second_vac = length(vec[which((kk-week_first_vec)>=dose_gap & vec==1)]) # check if any in cohort are due for second dose
    if(second_vac>0 && total_reserve_vac>0  && (total_reserve_vac)<=second_vac){  #if there are people waiting for second dose & there are doses available (exclude new doses set aside in current week)
      vec[which((kk-week_first_vec)>=dose_gap & vec==1)][1:total_reserve_vac] <- 2 #allocate available doses (in order of first dose receipt)
      week_second_vec[which((kk-week_first_vec)>=dose_gap)][1:total_reserve_vac] <- kk 
      total_reserve_vac <- total_reserve_vac - min(total_reserve_vac, second_vac)  #update reserve vaccine doses
    } else if (second_vac>0 && total_reserve_vac>0 && total_reserve_vac>second_vac) { #if there are more vaccines than people needing second dose 
      vec[which((kk-week_first_vec)>=dose_gap & vec==1)][1:second_vac] <- 2 #allocate available doses (in order of first dose receipt)
      week_second_vec[which((kk-week_first_vec)>=dose_gap)][1:second_vac] <- kk 
      total_reserve_vac <- total_reserve_vac - min(total_reserve_vac, second_vac) 
    } else {
      total_reserve_vac = total_reserve_vac 
    }
    
    #update matrices and vectors
    new_vac_vec[kk+1] <- new_vac
    new_reserve_vec[kk+1] <- new_reserve_vac
    new_distributed_vec[kk+1] <- new_distributed_vac
    total_reserve_vec[kk+1] <- total_reserve_vac
    vac_pop_mat[kk+1, ] <- vec
    
  }
  
  # account for waning protection if 2nd dose is delayed and delay from receipt of 1st dose to development of immunity
  tmp <- as.data.frame(vac_pop_mat) 
  tmp %>%
    mutate(across(, ~if_else(.==2, 999, .))) %>% 
    mutate(across(, ~if_else(.==1, cumsum(.),.)), 
           across(, ~if_else(.==0, 0, 
                             if_else(.==999, ve,
                                     if_else(.<=dose_gap & .<=delay & .>=1, 0,
                                             if_else(.<=(dose_gap + delay) & .>=1, ve*rel_protect_1_dose, ve*rel_protect_1_dose*waning_scalar^(.-dose_gap))))))) %>% 
    rowSums() -> rel_protect
  
  # waning protection and delay to start of protection, as above, but also protection drops to 0 if >6 weeks after 1st dose
  tmp %>%
    mutate(across(, ~if_else(.==2, -999, .))) %>% 
    mutate(across(, ~if_else(.==1, cumsum(.),.)), 
           across(, ~if_else((.==-999 & max(.)>dose_gap+3), 0, .)), # set protection to 0 for those who pass 6 weeks without receiving vaccine
           across(, ~if_else(.>0 & .>dose_gap+3, 0, .)), 
           across(, ~if_else(.==0, 0, 
                          if_else(.==-999, ve,
                                if_else(.<=dose_gap  & .<=delay & .>=1, 0,
                                    if_else(.<=(dose_gap + delay) & .>=1, ve*rel_protect_1_dose, ve*rel_protect_1_dose*waning_scalar^(.-dose_gap))))))) %>% 
    rowSums() -> rel_protect_drop
  
  # calculate gap between 2 doses
  tmp1 <- apply(vac_pop_mat==1, 2, which.max) 
  tmp2 <- apply(vac_pop_mat==2, 2, which.max) 
  tmp3 <- tmp2[tmp1>1 & tmp1<=(weeks_tot+1-dose_gap)] - tmp1[tmp1>1 & tmp1<=(weeks_tot+1-dose_gap)] #only calculate for those who received at least one dose and would have been eligible for 2nd dose by end of sim

  prop_second_dose <- length(tmp3[tmp3>1])/length(tmp3)
  prop_on_time <- length(tmp3[tmp3==dose_gap])/length(tmp3[tmp3>1])
  num_elig_second_dose <- length(tmp3)
  num_receive_second_dose <- length(tmp3[tmp3>1])
  num_second_dose_on_time <- length(tmp3[tmp3==dose_gap])
  num_any_dose <- length(tmp1[tmp1>1])
  num_within_6_weeks <- length(tmp3[tmp3<=dose_gap+3 & tmp3>1])

  #summarize outputs
  data.frame(vac_pop_mat) %>%  
    mutate(week = row_number()-1) %>%  
    pivot_longer(cols=-week, names_to="index", values_to="status") %>% 
    group_by(week, status) %>%
    summarize(tot = n(), .groups = 'drop') %>%
    pivot_wider(names_from="status", values_from=tot) %>%
    rename(unvac = `0`, 
           one_dose = `1`, 
           two_dose = `2`) %>%
    replace(is.na(.), 0) %>%
    mutate(tot_protection = ve*(one_dose * rel_protect_1_dose + two_dose)) %>% 
    bind_cols(new_vac = new_vac_vec, 
              new_distributed_vac = new_distributed_vec, 
              total_reserve_vac = total_reserve_vec, 
              new_reserve_vac = new_reserve_vec, 
              waning_protection = rel_protect, 
              waning_protection_drop = rel_protect_drop, 
              rr_inf_no_vac = rr_curve*params$base_inc) %>% 
    mutate(rr_inf_with_vac = (1 - waning_protection/pop_size) * rr_inf_no_vac , 
           risk_among_vac_base = (pop_size - unvac)*rr_inf_no_vac, 
           risk_among_vac = (pop_size - unvac - waning_protection) * rr_inf_no_vac, 
           vac_benefit = waning_protection * rr_inf_no_vac, 
           vac_benefit_drop = waning_protection_drop * rr_inf_no_vac, 
           prop_second_dose = round(prop_second_dose,2), 
           prop_on_time = round(prop_on_time,2), 
           num_any_dose = (params$pop_scalar/params$pop_size)*num_any_dose, 
           num_elig_second_dose = (params$pop_scalar/params$pop_size)*num_elig_second_dose, 
           num_receive_second_dose = (params$pop_scalar/params$pop_size)*num_receive_second_dose, 
           num_second_dose_on_time = (params$pop_scalar/params$pop_size)*num_second_dose_on_time, 
           num_within_6_weeks = (params$pop_scalar/params$pop_size)*num_within_6_weeks) -> output

  return(output)
  
}

plot_fun <- function(x){
  
  x%>%
    select(week, one_dose, two_dose, waning_protection, sim) %>%
    pivot_longer(cols=-c(week, waning_protection, sim), values_to="value", names_to="metric") %>%
    mutate(metric_num = ifelse(metric=="one_dose", 0.9, 1),
           metric = factor(metric, levels=c("one_dose", "two_dose"), labels=c("One dose", "Two doses"))
    ) %>%
    filter(week>0) %>%
    left_join(sims, by="sim") %>% 
    mutate(sim_name = factor(allocation_scheme, levels=c("fixed", "switch"), labels=c("Fixed", "Flexible")), 
           interact = paste0(sim_name, ", ", metric), 
           x_lab=1) %>% 
    ggplot(aes(x=interaction(x_lab, sim_name), y=(params$pop_scalar/params$pop_size)*value)) +
    geom_bar(aes(fill=interact, alpha=interact), stat="identity",   color="black") +
    geom_point(aes(y=(params$pop_scalar/params$pop_size)*waning_protection,color="Total effective population protection"),stroke=0.75, size=2, pch=21, fill="white") +
    facet_grid(~week, switch="x") +
    theme_classic() +
    theme(text = element_text(size=14), 
          legend.justification = "left",  
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          strip.background = element_blank()) +
    scale_color_manual(values="black") +
    scale_fill_manual(name = "Strategy and doses received", 
                      labels=c("Fixed, 1 dose", "Fixed, 2 doses", "Flexible, 1 dose", "Flexible, 2 doses"),
                      values=c("#EF8A62", "#EF8A62","#67A9CF", "#67A9CF")) +
    scale_alpha_manual(name = "Strategy and doses received", 
                       labels=c("Fixed, 1 dose", "Fixed, 2 doses", "Flexible, 1 dose", "Flexible, 2 doses"),
                       values = c(0.5, 1, 0.5, 1)) +
    labs(x="\nWeek", y="Number vaccinated (millions)\n",  color="")
  
}

rr_plot_fun <- function(x){

  x%>%
    select(week, sim, rr_inf_no_vac, rr_inf_with_vac) %>%
    mutate(pct_averted = 100*(rr_inf_no_vac-rr_inf_with_vac)/rr_inf_no_vac) %>%
    pivot_longer(cols=-c(week, sim), values_to="value", names_to="metric") -> x1
  
  x1 %>% 
    filter(week>0, metric=="rr_inf_no_vac") %>%
    left_join(sims, by="sim") %>%
    mutate(sim_name = allocation_scheme, 
           max_doses_in = paste0("Max weekly pop\nvaccinated: ", 100*max_doses_in, "%"), 
           inc_curve = inc )-> x2
    
  x1 %>% 
    filter(week>0, 
           metric=="rr_inf_with_vac") %>%
    left_join(sims, by="sim") %>%
    mutate(sim_name = allocation_scheme, 
           max_doses_in = paste0("Max weekly pop\nvaccinated: ", 100*max_doses_in, "%"), 
           inc_curve = inc ) %>% 
    ggplot(aes(x=week, y=value, color=sim_name)) +
    geom_line(lwd=1) + 
    geom_line(dat=x2, aes(y=value, lty="No vaccination"), color="black",  lwd=1) +
    facet_wrap(inc~.) +
    theme_bw() +
    scale_color_brewer(type="qual", palette=2) +
    scale_linetype_manual(values=c(2)) + 
    labs(x="\nWeek", y="Risk of infection\n", color="Allocation scheme", lty="") 
  
}

pct_averted_fun <- function(x){
  x %>%
    mutate(pct_averted = 100*(rr_inf_no_vac-rr_inf_with_vac)/rr_inf_no_vac) %>% 
    select(week, sim, pct_averted) %>% 
    left_join(sims, by="sim") %>%
    filter(week>0, inc==1) %>%
    mutate(sim_name = factor(allocation_scheme, levels=c("fixed", "switch"), labels=c("Fixed", "Flexible")), 
           scenario = scenario)
}


pct_averted_plot_fun <- function(x){
  x%>%
    mutate(pct_averted = 100*(rr_inf_no_vac-rr_inf_with_vac)/rr_inf_no_vac) %>% 
    select(week, sim, pct_averted) %>% 
    left_join(sims, by="sim") %>%
    filter(week>0, inc==1) %>%
    mutate(sim_name = factor(allocation_scheme, levels=c("fixed", "switch"), labels=c("Fixed", "Flexible"))) %>%
    ggplot(aes(x=week, y=pct_averted, color=sim_name)) +
      geom_point(size=2) +
      theme_classic() +
      theme(text = element_text(size=14)) +
      scale_color_manual(values=c("grey70", "grey30")) +
      labs(x="\nWeek", y="% of weekly symptomatic infections averted\nwith vaccination\n", color="Strategy")
  
}

delay_fun <- function(x, x1){
  x %>%
    group_by(sim) %>%
    filter(week==params$weeks_tot) %>%
    mutate(waning_protection = round((params$pop_scalar/params$pop_size) *  waning_protection, 2)) %>%
    select(sim, waning_protection) -> x2
  
  x%>%
    filter(week==0) %>%
    left_join(sims, by="sim") %>% 
    select("Allocation scheme" = allocation_scheme,
            inc = inc, 
           "# receiving at least one dose (millions)" = num_any_dose,
           "# eligible to receive second dose (millions)" = num_elig_second_dose,
           "# receiving second dose (millions)" = num_receive_second_dose, 
           "# receiving second dose within 6 weeks (millions)" = num_within_6_weeks, 
           "# receiving second dose on week 3 (millions)" = num_second_dose_on_time, 
           sim) %>%
    left_join(x2, by="sim") %>%
    left_join(x1, by=c("Allocation scheme"="scheme", "inc")) %>% 
    mutate(inc = factor(inc, levels=seq(1:3), labels=c("Stable", "Increasing", "Peaking")), 
           scenario = scenario, 
           `Allocation scheme` = if_else(`Allocation scheme`=="switch", "Flexible", "Fixed")) %>%
    select(scenario, everything(), -max_doses_in, -abs_diff, -sim) %>%
    rename(Scenario = scenario, 
           "Incidence trend" = inc,
           "% averted (flexible compared to fixed)" = pct_diff, 
           "Total effective vaccine protection (millions)" =waning_protection) %>%
    arrange(`Incidence trend`)
    
}

delay_fun_drop <- function(x, x1){
  x %>%
    group_by(sim) %>%
    filter(week==params$weeks_tot) %>%
    mutate(waning_protection = round((params$pop_scalar/params$pop_size) * waning_protection_drop, 2)) %>%
    select(sim, waning_protection) -> x2
  x%>%
    filter(week==0) %>%
    left_join(sims, by="sim") %>% 
    select("Allocation scheme" = allocation_scheme,
           inc = inc, 
           "# receiving at least one dose (millions)" = num_any_dose,
           "# eligible to receive second dose (millions)" = num_elig_second_dose,
           "# receiving second dose (millions)" = num_receive_second_dose, 
           "# receiving second dose within 6 weeks (millions)" = num_within_6_weeks, 
           "# receiving second dose on week 3 (millions)" = num_second_dose_on_time, 
           sim) %>%
    left_join(x2, by="sim") %>%
    left_join(x1, by=c("Allocation scheme"="scheme", "inc")) %>% 
    mutate(inc = factor(inc, levels=seq(1:3), labels=c("Stable", "Increasing", "Peaking")), 
           scenario = scenario, 
           `Allocation scheme` = if_else(`Allocation scheme`=="switch", "Flexible", "Fixed")) %>%
    select(scenario, everything(), -max_doses_in, -abs_diff, -sim) %>%
    rename(Scenario = scenario, 
           "Incidence trend" = inc,
           "% averted (flexible compared to fixed)" = pct_diff, 
           "Total effective vaccine protection (millions)" =waning_protection) %>%
    arrange(`Incidence trend`)
  
}

cum_benefit_fun <- function(x){
  x %>% 
    group_by(sim) %>%
    summarize(cum_inf_risk = sum(risk_among_vac_base), 
              cum_benefit_vac = sum(vac_benefit), .groups="drop") %>%
    left_join(sims, by="sim") %>%
    select(cum_benefit_vac, allocation_scheme, max_doses_in, inc) %>%
    pivot_wider(id_cols=c(max_doses_in, inc), names_from=allocation_scheme, values_from=cum_benefit_vac) %>%
    mutate(pct_diff = round(100*(switch-fixed)/fixed,1), 
           abs_diff = round((switch-fixed),1)) %>%
    select(inc, max_doses_in, pct_diff, abs_diff) %>%
    mutate(scheme="switch")
  
}

cum_benefit_drop_fun <- function(x){
  x %>% 
    group_by(sim) %>%
    summarize(cum_benefit_vac = sum(vac_benefit_drop), .groups="drop") %>%
    left_join(sims, by="sim") %>%
    select(cum_benefit_vac, allocation_scheme, max_doses_in, inc) %>%
    pivot_wider(id_cols=c(max_doses_in, inc), names_from=allocation_scheme, values_from=cum_benefit_vac) %>%
    mutate(pct_diff = round(100*(switch-fixed)/fixed,1), 
           abs_diff = round((switch-fixed),1)) %>%
    select(inc, max_doses_in, pct_diff, abs_diff) %>%
    mutate(scheme="switch")
  
}
