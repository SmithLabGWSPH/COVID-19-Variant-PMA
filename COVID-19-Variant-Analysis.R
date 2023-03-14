library(dplyr)
library(tidyr)
library(lubridate)
library(haven)
library(readxl)
library(ggplot2)
library(stringr)


# Create output folders

dir.create("data_out", showWarnings = FALSE) # Creating a folder data_out if it doesn't exist
dir.create("plots", showWarnings = FALSE) # Creating a folder plots if it doesn't exist
dir.create("plots/variants_nextstrain", showWarnings = FALSE)

###########################
# Prepare Nextstrain data #
###########################
# Read in nextstrain raw data

nextstrain_frequency <- read.csv('data/nextstrain_frequencies_2022oct05.csv')
nextstrain_metadata <- read.csv('data/nextstrain_metadata_2022oct05.csv')
nextstrain_pivots <- read.csv('data/nextstrain_pivotdates_2022oct05.csv', header = FALSE)

# Join/merge nextstrain data sets to create strain frequency by week and by country

nextstrain_df <- nextstrain_metadata %>% select(strain, date, country, clade_membership)
nextstrain_df <- left_join(nextstrain_frequency, nextstrain_df, by="strain")

pivot_df <- data.frame(week = 0:(length(nextstrain_pivots$V1)-1), pivot = nextstrain_pivots$V1)
#converting the dates in decimal in nextstrain dataset to readable dates.

pivot_df <- pivot_df %>% mutate(date = as.Date(date_decimal(pivot))) 
nextstrain_df <- left_join(nextstrain_df, pivot_df, by="week")

# Aggregate clades into variants

nextstrain_df <- nextstrain_df  %>%
  mutate(variant = case_when(clade_membership %in% 
                               c("19A", "19B", "20A", "20B", "20C", "20D", "20F", "20G", "20E (EU1)") ~ "Pre_alpha",
                             clade_membership %in% 
                               c("20I (Alpha, V1)") ~ "Alpha",
                             clade_membership %in%
                               c("20H (Beta, V2)") ~ "Beta", 
                             clade_membership %in%
                               c("21B (Kappa)") ~ "Kappa", 
                             clade_membership %in%
                               c("21D (Eta)") ~ "Eta",
                             clade_membership %in% 
                               c("21E (Theta)") ~ "Theta",
                             clade_membership %in% 
                               c("20J (Gamma, V3)") ~ "Gamma",
                             clade_membership %in% 
                               c("21A (Delta)", "21I (Delta)", "21J (Delta)") ~ "Delta",
                             clade_membership %in% 
                               c("21C (Epsilon)") ~ "Epsilon",
                             clade_membership %in% 
                               c("21F (Iota)") ~ "Iota",
                             clade_membership %in% 
                               c("21G (Lambda)") ~ "Lambda",
                             clade_membership %in% 
                               c("21H (Mu)") ~ "Mu",
                             clade_membership %in% 
                               c("21K (Omicron)", "21L (Omicron)", 
                                 "21M (Omicron)", "22A (Omicron)","22B (Omicron)", 
                                 "22C (Omicron)", "22D (Omicron)") ~ "Omicron"))


# Sum frequencies of individual clades into frequencies by variant (`freq`)

nextstrain_df  <- nextstrain_df  %>%
  group_by(country, week, variant) %>%
  summarize(freq = sum(frequency))

# Compute relative/normalized frequency (`nor`) per variant per week per country, from absolute frequency (prevalence)

nextstrain_df  <- nextstrain_df  %>%
  group_by(country, week) %>% 
  mutate(nor = (freq/sum(freq))) %>%
  filter(country %in% c("Kenya", "Spain", "Hong Kong", "Italy",
                        "USA", "Mexico", "Switzerland", "France")) 

# Select only the last week of each month

nextstrain_df <- left_join(nextstrain_df, pivot_df, by="week")
nextstrain_df <- nextstrain_df %>%
  mutate(date_ym = format(date, format= "%y-%m"))
nextstrain_df <- nextstrain_df %>%
  group_by(country, date_ym) %>%
  mutate(max_date = max(date)) %>% 
  filter(date == max_date)

# Read in individual patient data (IPD)

compiled_ipd <- read_dta('data/Compiled_IPD_Data_221130.dta')
site_metadata <- read.csv('data/country_id_ipd_2.csv') %>%
  mutate(site_id = as.integer(site_id))

####
# Read in Aggregate data from Switzerland and France
#
# Switzerland and France data are in multi-sheet Excel workbooks, with the table for each outcome on a separate sheet.
# The name of each sheet is the outcome for that sheet
####

swiss_sheets_vec <- excel_sheets('data/Variant-Swiss-France-Data_20221208/TABLE_outcomes_switzerland.xlsx') # creates a vector.
swiss_df <- data.frame(month_year = {}, Events={},	Missing={},	Total={})
for (sheet in swiss_sheets_vec) {
  sheet_df <- read_excel('data/Variant-Swiss-France-Data_20221208/TABLE_outcomes_switzerland.xlsx', 
                         sheet=sheet, 
                         skip=2, 
                         col_names=c("month_year", "Events", "Missing", "Total"))
  sheet_df <- sheet_df %>% dplyr::mutate(Outcome=sub("swiss - ", "", sheet),
                                         country="Switzerland")
  swiss_df <- rbind(swiss_df, sheet_df)
}

# France
france_sheets_vec <- excel_sheets('data/Variant-Swiss-France-Data_20221208/TABLE_outcomes_france.xlsx') # creates a vector.
france_df <- data.frame(month_year = {}, Events={},	Missing={},	Total={})
for (sheet in france_sheets_vec) {
  sheet_df <- read_excel('data/Variant-Swiss-France-Data_20221208/TABLE_outcomes_france.xlsx', 
                         sheet=sheet, 
                         skip=2, 
                         col_names=c("month_year", "Events", "Missing", "Total"))
  sheet_df <- sheet_df %>% dplyr::mutate(Outcome=sub("france - ", "", sheet),
                                         country="France")
  france_df <- rbind(france_df, sheet_df)
}

# Concatenate Switzerland and France data

swiss_france_df <- rbind(swiss_df, france_df)

# Transform dates to match format of IPD 

swiss_france_df <- swiss_france_df %>%
  mutate(month_year = lubridate::parse_date_time(month_year, 'my'),
         covid_date_ym = format(month_year, format= "%y-%m")) %>% 
  select(-month_year) %>%
  rename(outcome=Outcome,
         events1=Events,
         eventsm=Missing) %>%
  mutate(site_name = country,
         site_id=16,
         author_year = "Favre, Panchaud, 2021")

# Add 2 rows into site_metadata for France and Switzerland
# 
site_metadata <- site_metadata %>%
  rbind(list(16, "", "Switzerland", "Switzerland_COVIPREG", "Favre, Panchaud, 2022")) %>%
  rbind(list(16, "", "France", "France_COVIPREG", "Favre, Panchaud, 2022")) %>%
  filter(site_name != "COVIPREG")

# REMOVING PAKISTAN (AMHANII), BRAZIL (REBRACO) FROM THE DATA:THE OVERLAP PAST PreALPHA for rebraco IS NOT ENOUGH.

compiled_ipd <- compiled_ipd %>% filter(site_id !=38)
compiled_ipd <- compiled_ipd %>% filter(site_id !=36)

# Deriving just the year and month because it will be combined with 
# Nextstrain data, which we are using at a monthly level.
# This is also why we drop IPD that lacks a covid date
compiled_ipd <- compiled_ipd %>%
  mutate(covid_date_ym = format(date_onset_or_test, format="%y-%m")) %>%
  filter(!is.na(covid_date_ym))

compiled_ipd <- left_join(compiled_ipd, site_metadata, by="site_id")

#######
# The next section is deriving which months of the pandemic
# each study's data covers, in order to draw the plots showing
# which months each study covers, plotted against the variants in that country per month
######
study_months <- compiled_ipd %>%
  select(site_id, city, country, site_name, author_year, date_onset_or_test, covid_date_ym) %>%
  group_by(country, author_year, covid_date_ym) %>%
  summarize(Events = n())

## add in Switzerland and France
study_months_sf <- swiss_france_df %>%
  group_by(country, covid_date_ym, author_year) %>% 
  summarize(Events = sum(events1))

study_months <-
  rbind(study_months, study_months_sf) %>%
  mutate(covid_date_label = paste0("20", covid_date_ym, "-01")) %>%
  mutate(pandemic_month = interval("2019-12-01", covid_date_label) %/% months(1)) %>%  # "modulo" - give number of whole months
  select(-covid_date_ym) %>%
  mutate(study_type = "RFA",
         study_type_full = "Risk Factor Analysis")

#### CONTINUE WORKING HERE
#### "df8" is nextstrain_df, which has been created
#### the other DF we need is studies_df
#### Proceed borrowing from line 510 "Nextstrain figure"

study_months_agg_df <- study_months %>%
  filter(Events >= 1) %>%                          # only include months that had >=1 events for that study
  group_by(country, author_year, study_type) %>%
  summarize(min_month = min(pandemic_month),       # earliest month in the study with data (we're ignoring gaps)
            max_month = max(pandemic_month)) %>%   # latest month in the study with data   (we're ignoring gaps)
  group_by(country) %>%
  mutate(country_study_number = 1:n(),  # number the studies from the same country 1..n
         country_total_studies = n())   # total number of studies in the country


# vertical gap between bars representing studies
# (if there are 2 or more studies in the country)
gap_betw_studies <- 0.15

studies_variants_plot <-
  ggplot() +
  geom_bar(data = nextstrain_df, stat="identity",
           mapping = aes(x = date_ym, fill = variant, y=nor)) +
  scale_fill_manual(values = c("Alpha" = "deepskyblue2", "Pre_alpha" = "lightsteelblue1", "Beta" = "darkseagreen1",
                               "Delta" = "gold", "Epsilon" = "burlywood", "Eta" = "cadetblue2", "Gamma" = "darkgoldenrod1",
                               "Iota" = "cyan3", "Lambda" = "gray70", "Mu" = "peachpuff1", "Omicron" = "plum2",
                               "Theta" = "tomato")) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
        axis.text.y = element_text(size = 6)) + 
  xlab("Year-Month") + 
  ylab("Relative Frequency") +
  geom_rect(data = study_months_agg_df,
            mapping=aes(xmin=min_month-0.5,
                        xmax=max_month+0.5,
                        ymin= 0.5 - 0.5*gap_betw_studies*(country_total_studies-1) + gap_betw_studies*(country_study_number - 1) - 0.005,
                        ymax= 0.5 - 0.5*gap_betw_studies*(country_total_studies-1) + gap_betw_studies*(country_study_number - 1) + 0.005)) +
  geom_text(data = study_months_agg_df,
            mapping=aes(x = 0.5*(min_month + max_month),
                        y = 0.5 - 0.5*gap_betw_studies*(country_total_studies-1) + gap_betw_studies*(country_study_number - 1) + 0.05,
                        label = author_year),
            #   label = ifelse(is.na(site_name),
            #                 study_type,
            #                paste(site_name, study_type))),
            size = 2) +
  #guides(fill = guide_legend(override.aes = list(labels = c("A", "B")))) +
  facet_wrap(~country, nrow = 30, ncol = 1, scale = "free_x")

today_date <- format(Sys.time(), "%Y%m%d")
ggsave(plot = studies_variants_plot,
       filename = paste0('plots/variants_nextstrain/covid_variant_country_',
                         today_date,'.pdf'),
       width = 6, height = 20, units = "in")
  
##################
# STUDY ANALYSIS #
##################

# Remove twins and triplets

compiled_ipd <- compiled_ipd %>%
  dplyr::filter(twin2!=1) %>%
  dplyr::filter(tripletnum==0)

compiled_ipd <- compiled_ipd %>%
  select(site_id, city, country, site_name, author_year, person_id,
         covid_date_ym,
         twin, triplet, twin2, tripletnum, 
         mat_age, matage_cat, matage19, matage2024, matage2529, matage3034, 
         matage3539, matage40, matagem, matage1519, matage4045,matage35_older,
         bmi_pre, bmi_pre_cat, bmi_prepreg_under, bmi_prepreg_normal, 
         bmi_prepreg_over, bmi_prepreg_obese, bmi_prepreg_miss,
         # maternal outcomes: 
         diab_prepreg1, diab_prepreg2, diab_prepregm, 
         hyperten_prepreg1, hyperten_prepreg2,hyperten_prepregm, 
         gestage_covid_tri1, gestage_covid_tri2, gestage_covid_tri3,
         gestage_covid_postp, gestage_covid_trim, date_onset_or_test,
         covid_sympdate, covid_pos, covid_hosp1, covid_hosp0, covid_hospm,
         fetuses, parity, gravidity, weight_pre,
         prdeath1, prdeath0, prdeathm,
         ICUadmit1, ICUadmit0, ICUadmitm,
         critcare1,critcare0, critcarem, 
         ventilation1, ventilation0,ventilationm, 
         pneumonia1, pneumonia0, pneumoniam,
         prdeath1, prdeath0, prdeathm, 
         place_abrupt1, place_abrupt0,place_abruptm, 
         preterm_labor1, preterm_labor0, preterm_laborm,
         haemorrhage1, haemorrhage0, haemorrhagem, 
         embolicdz1, embolicdz0, embolicdzm,
         preeclampsia1, preeclampsia0, preeclampsiam, 
         eclampsia1, eclampsia0, eclampsiam, 
         pre_or_eclampsia1, pre_or_eclampsia0, pre_or_eclampsiam, 
         hpd_any1, hpd_any0, hpd_anym, 
         hpd_postcovid1, hpd_postcovid0, hpd_postcovidm, 
         c_section1, c_section0, c_sectionm, 
         c_intrap1, c_intrap0, c_intrapm,
         
         # fetal outcomes: 
         stillbirth_site1, stillbirth_site0, stillbirth_sitem,
         stillbirth_281, stillbirth_280, stillbirth_28m,
         perinatal_death281, perinatal_death280, perinatal_death28m,
         perinatal_deathsite1, perinatal_deathsite0, perinatal_deathsitem, 
         nicu1, nicu0, nicum, 
         sga1, sga0, sgam, 
         extremesga1, extremesga0, extremesgam,
         lowbirthweight1, lowbirthweight0, lowbirthweightm,
         verylowbirthweight1, verylowbirthweight0, verylowbirthweightm,
         sab_site1, sab_site0, sab_sitem,
         tab_site1, tab_site0, tab_sitem,
         preterm1, preterm0, pretermm,
         verypreterm1, verypreterm0, verypretermm,
         neonatal_death1, neonatal_death0, neonatal_deathm,
         earlyneo_death1, earlyneo_death0, earlyneo_deathm,
         
         # denominators:
         denom_allpreg,denom_endpreg,denom_endpreg37,denom_totalbirth28wk,
         denom_livebirth,denom_livebirth34, denom_livebirth37) 

####
# Look at these results, determine what to do, then ultimately remove this code
###

check_df <- compiled_ipd %>%
  select(preterm1, verypreterm1, preterm_labor1,
         denom_livebirth37, denom_livebirth34, denom_endpreg37) %>%
  mutate(preterm21 = if_else((preterm1==1 & denom_livebirth37==1), 1, 0),
         verypreterm21 = if_else((verypreterm1==1 & denom_livebirth34==1), 1, 0),
         preterm_labor21 = if_else((preterm_labor1==1 & denom_endpreg37==1), 1,0))

check_preterm_tab <-
  table(preterm1 = compiled_ipd$preterm1, denom_livebirth37 = compiled_ipd$denom_livebirth37)
check_preterm_tab

check_verypreterm_tab <-
  table(verypreterm1 = compiled_ipd$verypreterm1, denom_livebirth34 = compiled_ipd$denom_livebirth34)
check_verypreterm_tab

check_preterm_labor_tab <-
  table(preterm_labor1 = compiled_ipd$preterm_labor1, denom_endpreg37 = compiled_ipd$denom_endpreg37)
check_preterm_labor_tab

###### Resume here

table1_df <- compiled_ipd %>%
  group_by(city, country) %>%
  summarize(mean_age  = mean(mat_age, na.rm = TRUE), 
            sd_age = sd(mat_age, na.rm = TRUE),
            bmi_mean = mean(bmi_pre, na.rm = TRUE),
            bmi_under = sum(bmi_prepreg_under),
            bmi_normal = sum(bmi_prepreg_normal),
            bmi_over = sum(bmi_prepreg_over),
            bmi_obese = sum(bmi_prepreg_obese),
            bmi_missing = sum(bmi_prepreg_miss),
            tot_preg = sum(denom_allpreg),
            livebirths = sum(denom_livebirth),
            ga_tri1 = round((sum(gestage_covid_tri1)/tot_preg)*100, 0),
            ga_tri2 = round((sum(gestage_covid_tri2)/tot_preg)*100, 0),
            ga_tri3 = round((sum(gestage_covid_tri3)/tot_preg)*100, 0),
            ga_pp = round((sum(gestage_covid_postp)/tot_preg)*100, 0),
            ga_unknown = round((sum(gestage_covid_trim)/tot_preg)*100, 0),
            hospitalized = round(sum(covid_hosp1)/tot_preg*100, 0),
            icu_admit = round(sum(ICUadmit1)/tot_preg*100, 0),
            twins = sum(twin2),
            triplets = sum(tripletnum))

# NOTE: THIS COMPOSITE OUTCOME CAN'T BE DONE ON AGGREGATE DATA (SWISS FRANCE) B/C EXAMPLE: ICU ADMISSION EVENTS =3, CRITICAL CARE EVENTS = 3, WE DON'T KNOW IF THOSE TWO OUTCOMES HAPPENED TO THE SAME 3 WOMEN OR 6 DIFFERENT WOMEN. 
# So the composite outcome will be missing for Switzerland and France.
#
# Note on the code below:
# If ANY of the outcome variables (for example, place_abrupt1, etc.) are 1, then the composite outcome will be 1
# EVEN if other variables are NA  (because NA | NA |...| 1 == TRUE).  This is what we want.

compiled_ipd <- compiled_ipd %>% 
  mutate(mat_composite_outcome1 = if_else((place_abrupt1==1 |
                                             preterm_labor1==1 | 
                                             haemorrhage1==1 | 
                                             embolicdz1==1 | 
                                             preeclampsia1==1 | eclampsia1==1 | pre_or_eclampsia1==1 |
                                             hpd_any1==1 | 
                                             c_section1==1 | 
                                             c_intrap1==1), 1, 0), 
         mat_composite_outcome0 = 1-mat_composite_outcome1, 
         mat_composite_outcomem = 0)  %>%
  mutate(fetal_composite_outcome1 = if_else((stillbirth_281==1 | 
                                               neonatal_death1==1 | earlyneo_death1==1 | 
                                               perinatal_death281==1 | 
                                               nicu1==1 |
                                               lowbirthweight1==1 | verylowbirthweight1==1 | 
                                               preterm1==1 |
                                               verypreterm1==1 | 
                                               sga1==1 | extremesga1==1), 1, 0), 
         fetal_composite_outcome0 = 1-fetal_composite_outcome1, 
         fetal_composite_outcomem = 0)

###################################################################
# Calculating dominant strain in per month, per country           #
# Defined as the strain whose frequency is >= the given threshold #
###################################################################

dom_threshold <- 0.50
nextstrain_agg_df <- nextstrain_df %>%
  group_by(country, date_ym) %>%
  arrange(desc(freq)) %>% # we first ordered the frequency and the highest is on top. 
  mutate(dominant_variant = if_else(first(freq)==0,  # if the highest non-missing frequency strain is 0
                                    "no_strain",     # then set dominant variant as "no_strain"
                                    # Otherwise, if we have some numbers: 
                                    ifelse(!is.na(first(nor)) & first(nor)>=dom_threshold, # if we DO have non-zero frequencies
                                           first(variant), # and the top frequency is > dom_threshold, then this is our dominant variant 
                                           "mixed"))) %>%  # otherwise set dominant variant to "mixed>
  ungroup() %>%
  arrange(country, date_ym) # arrange the strains by country and date order. 

nextstrain_agg_wide_df <- nextstrain_agg_df %>% 
  select(-c(freq, week, date, max_date, pivot)) %>% 
  pivot_wider(names_from = variant,
              values_from = nor,
              values_fill = 0,
              names_prefix = "freq_")

compiled_ipd <- left_join(compiled_ipd, nextstrain_agg_wide_df,
                          by = join_by('country' == 'country',
                                       'covid_date_ym' == 'date_ym'))

outcome_events <- compiled_ipd %>%
  pivot_longer(cols=c(prdeath1, prdeath0, prdeathm,
                      covid_hosp1, covid_hosp0, covid_hospm,
                      ICUadmit1, ICUadmit0, ICUadmitm,
                      critcare1,critcare0, critcarem, 
                      ventilation1, ventilation0,ventilationm, 
                      pneumonia1, pneumonia0, pneumoniam,
                      prdeath1, prdeath0, prdeathm, 
                      place_abrupt1, place_abrupt0,place_abruptm, 
                      preterm_labor1, preterm_labor0, preterm_laborm,
                      haemorrhage1, haemorrhage0, haemorrhagem, 
                      embolicdz1, embolicdz0, embolicdzm,
                      preeclampsia1, preeclampsia0, preeclampsiam, 
                      eclampsia1, eclampsia0, eclampsiam, 
                      pre_or_eclampsia1, pre_or_eclampsia0, pre_or_eclampsiam, 
                      hpd_any1, hpd_any0, hpd_anym, 
                      hpd_postcovid1, hpd_postcovid0, hpd_postcovidm, 
                      c_section1, c_section0, c_sectionm, 
                      c_intrap1, c_intrap0, c_intrapm,
                      # fetal outcomes: 
                      # stillbirth_site1, stillbirth_site0, stillbirth_sitem,
                      stillbirth_281, stillbirth_280, stillbirth_28m,
                      perinatal_death281, perinatal_death280, perinatal_death28m,
                      
                      # perinatal_deathsite1, perinatal_deathsite0, perinatal_deathsitem, 
                      nicu1, nicu0, nicum, 
                      sga1, sga0, sgam, 
                      extremesga1, extremesga0, extremesgam,
                      lowbirthweight1, lowbirthweight0, lowbirthweightm,
                      verylowbirthweight1, verylowbirthweight0, verylowbirthweightm,
                      #  sab_site1, sab_site0, sab_sitem,
                      #   tab_site1, tab_site0, tab_sitem,
                      preterm1, preterm0, pretermm, 
                      verypreterm1, verypreterm0, verypretermm, 
                      neonatal_death1, neonatal_death0, neonatal_deathm,
                      earlyneo_death1, earlyneo_death0, earlyneo_deathm,
                      mat_composite_outcome1, mat_composite_outcome0, mat_composite_outcomem,
                      fetal_composite_outcome1, fetal_composite_outcome0, fetal_composite_outcomem),
               
               # denominators:
               #denom_allpreg,denom_endpreg,denom_endpreg37,denom_totalbirth28wk,
               # denom_livebirth,denom_livebirth34, denom_livebirth37),
               names_to = "outcome", values_to = "events")

outcome_events <- outcome_events %>%
  group_by(site_name, author_year, country, dominant_variant, outcome, covid_date_ym) %>%
  summarize(events = sum(events, na.rm = FALSE))

# Mutate - remove the 0,1,m modifier and put in a new column.
outcome_events <- outcome_events %>%
  mutate(outcome_mod = str_sub(outcome, -1), # from the first position all to way to second from the end.
         outcome = str_sub(outcome, 1,-2)) # the last char. 

outcome_events <- outcome_events %>%
  pivot_wider(names_from = outcome_mod, names_prefix = "events", values_from = events) %>%
  mutate(Total = events0 + events1 + eventsm) 

# Prepare the Switzerland/France data with the same
# event-count aggregation

swiss_france_df <- left_join(swiss_france_df,
                             nextstrain_agg_wide_df,
                             by=join_by('country' == 'country',
                                        'covid_date_ym' == 'date_ym'))

swiss_france_df <- swiss_france_df %>%
  select(country, site_name, author_year, outcome,
         covid_date_ym, dominant_variant,
         events1, eventsm, Total)

# Mutate - remove the 0,1,m modifier and put in a new column
swiss_france_df <- swiss_france_df %>%
  mutate(outcome_mod = str_sub(outcome, -1), # from the first position all to way to second from the end.
         outcome = str_sub(outcome, 1,-2)) # the last char. 

# Merge in Switzerland/France data with aggregated IPD data
outcome_events <- rbind(outcome_events, swiss_france_df)

outcome_events <- outcome_events %>%
  mutate(incidence = events1/Total)

# Remove outcome for a given country if >25% of data for
# that country and outcome is missing
missing_threshold <- 0.25

outcome_events <- outcome_events %>% group_by(country, outcome) %>% 
  mutate(outcome_missing = sum(eventsm),
         outcome_total = sum(Total),
         pct_missing = outcome_missing/outcome_total) %>%
  select(-outcome_missing, -outcome_total) %>% 
  filter(pct_missing < missing_threshold)

################
# Analysis
###############

# Functions to be used

### used this link to calc: https://www.medcalc.org/manual/relative-risk-odds-ratio.php
### function for calculating RR (not logRR)
rr_ff <- function(pos_test, pos_ref, neg_test, neg_ref, adj=0.5){
  useadj <- pos_test==0 | neg_test==0 | pos_ref==0 | neg_ref==0
  
  a <- pos_test+useadj*adj
  b <- neg_test+useadj*adj
  c <- pos_ref+useadj*adj
  d <- neg_ref+useadj*adj
  return((a/(a+b))/((c/(c+d))))
}

### this calculates SE for logRR not just RR - so this is logSE
se_logrr_ff <- function(pos_test, pos_ref, neg_test, neg_ref, adj=0.5){
  useadj <- pos_test==0 | neg_test==0 | pos_ref==0 | neg_ref==0
  a <- pos_test+useadj*adj
  b <- neg_test+useadj*adj
  c <- pos_ref+useadj*adj
  d <- neg_ref+useadj*adj
  return(sqrt((1/a)+(1/c)-(1/(a+b))-(1/(c+d))))
}

# Pick up at line 755