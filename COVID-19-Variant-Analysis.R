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


outcome_events <- outcome_events %>%
  mutate(events_neg = Total - events1 - eventsm)

# vector of all strain names
strain_vec <- unique(outcome_events$dominant_variant)
strain_vec <- strain_vec[!is.na(strain_vec) & strain_vec != "no_strain"]

# Keep track of statistics for each study/outcome combination
studies_outcomes <- data.frame()

# iterate through pairs of strains
for (ref.strain in strain_vec) {
  for (test.strain in strain_vec) {
    if (ref.strain == test.strain) {
      next # skip when the two strains are the same
    }
    study_rows <- outcome_events %>%
      filter(dominant_variant %in% c(ref.strain, test.strain)) %>%
      group_by(site_name, author_year, dominant_variant, outcome) %>%
      #  dplyr::filter(outcome %in% c("matdeath", "preterm_labor")) %>% 
      summarise(pos = sum(events1), 
                neg = sum(events_neg)) %>%
      mutate(dominant_variant = if_else(dominant_variant==ref.strain,
                                        "reference",
                                        "test")) %>%
      pivot_wider(names_from = dominant_variant, values_from = c("pos", "neg")) %>%
      drop_na() %>% # drops countries that don't have one or the other dominant strain.
      mutate(ref.strain=ref.strain, test.strain=test.strain) %>%
      mutate(logRR = if_else(pos_reference==0 & pos_test==0,
                             as.numeric(NA),
                             log(rr_ff(pos_ref=pos_reference, pos_test=pos_test, neg_ref=neg_reference, neg_test=neg_test))),
             logSE = if_else(pos_reference==0 & pos_test==0,
                             as.numeric(NA),
                             se_logrr_ff(pos_ref=pos_reference, pos_test=pos_test, neg_ref=neg_reference, neg_test=neg_test)))
    
    studies_outcomes <- rbind(studies_outcomes, study_rows)
  }
}


# Function to perform meta-analysis (and save/create forest plot)
# for a given pair of strains

metagen_fun <- function(refstr, teststr) {
  strainpair_df <- studies_outcomes %>%
    filter(ref.strain==refstr & test.strain==teststr)
  # relabel test and ref pos/neg columns to, for example, Omicron+, Omicron-, Delta+, Delta-
  colnames(strainpair_df)[colnames(strainpair_df) == "pos_reference"] <- paste0(refstr, "+")
  colnames(strainpair_df)[colnames(strainpair_df) == "neg_reference"] <- paste0(refstr, "-")
  colnames(strainpair_df)[colnames(strainpair_df) == "pos_test"] <- paste0(teststr, "+")
  colnames(strainpair_df)[colnames(strainpair_df) == "neg_test"] <- paste0(teststr, "-")
  
  if(nrow(strainpair_df)==0){ # when a certain pair of strains is not present, we wnat it to skip, e.g. Delta-Epsilon is missing when dominant strain defined as >=70% prevalence. 
    print("     No such strain pair found - skipping meta-analysis")
    return()
  }
  
  ### metagen 
  mymeta <- metagen(logRR, logSE,
                    studlab = paste0(site_name, " (", author_year, ")"), data=strainpair_df,
                    fixed = FALSE, random=TRUE,
                    subset = NULL, sm = "RR",
                    method.tau = "ML", # Default Restricted maximum likelihood estimator (REML) did not always converge so using ML. 
                    subgroup = strainpair_df$outcome,
                    title = "Variant")
  
  # added in the continuity correction by hand - coded it up above
  # manual process in STATA - replaced 0 with the inverse of the opposite arm. can try this code. Making the continuity correction smaller helps but not always. Pooled absolute risks and pooled relative risks that were illogical, absolute risk was high for one group but relative risk was lower for that group.   If there are 0 events for delta and alpha, what to do with that estimate. 
  # summary (mymeta)
  
  
  # Extracting the follwing statistics:
  #
  # TE.random.w (Estimated effect in subgroup (random effect model))
  # lower.random.w (CI for the random effects)
  # upper.random.w (upper CI for the random effects)
  # k.study.w has number of studies that are included for each random effects model for each outcome given a ref strain and test strain. 
  # for more info: ?metagen --> then click on "meta.object".
  # http://127.0.0.1:57661/help/library/meta/help/meta-object
  
  strainpair_meta_results <- data.frame(Outcome=mymeta$bylevs, Strain=teststr, Ref_Strain=refstr, 
                                        RR=exp(mymeta$TE.random.w), LCI=exp(mymeta$lower.random.w), UCI=exp(mymeta$upper.random.w), n_study=mymeta$k.study.w)
  
  # generate forest plot and save
  
  day_string <- format(Sys.time(), "%Y%m%d")
  dir_folder_name <- paste0(day_string,"_domvar_mixed_",as.character(dom_threshold*100))
  dir.create(paste0("plots/",dir_folder_name), recursive = TRUE)
  pdf(file=paste0("plots/",dir_folder_name,"/VariantStudy_MetaForestPlot_teststr_", teststr, "_refstr_", refstr, ".pdf"), width=10, height=65) # pdf saving has to go before making the plot and at the end have to say dev.off()
  
  forest(mymeta, sortvar=logRR,
         leftcols = c("studlab",
                      paste0(teststr, "+"), paste0(teststr, "-"),
                      paste0(refstr, "+"), paste0(refstr, "-")),
         lab.NA.effect = "NA",
         colgap.forest.left = "1mm",
         fontsize = 10,
         fs.hetstat = 9,
         header=TRUE,
         overall=FALSE,
         col.diamond = "#33CC66",
         col.study= "#333333")
  
  dev.off()
  return(strainpair_meta_results)
}

# Perform meta-analysis (and create/save forest plots) for each pair of strains

heatmap_data_df <- data.frame(Outcome = {}, Strain = {}, Ref_Strain = {}, RR = {}, LCI = {}, UCI = {}, n_study = {})

for (ref.strain in strain_vec) {
  for (test.strain in strain_vec) {
    if(ref.strain==test.strain) {
      next # skip this iteration when the two strains are the same. 
    }
    print(paste(ref.strain, test.strain))
    strainpair_results <- metagen_fun(ref.strain, test.strain)
    heatmap_data_df <- rbind(heatmap_data_df, strainpair_results)
  }
}

heatmap_data_df <- heatmap_data_df %>% 
  mutate(RR_CI_lab = if_else(is.na(RR),
                             "n=0",
                             paste0(format(round(RR, 2), nsmall = 2),
                                    " (",
                                    format(round(LCI, 2), nsmall = 2),
                                    ", ",
                                    format(round(UCI, 2), nsmall = 2),
                                    "), n=",
                                    n_study)))

est_table_df <- heatmap_data_df %>%
  select(-RR, -LCI, -UCI, -n_study) %>%
  pivot_wider(names_from = Strain, values_from = RR_CI_lab, names_prefix = "Strain_") %>% 
  pivot_longer(cols = starts_with("Strain_"), names_to = "Strain", values_to = "RR_CI_lab")  %>%
  mutate(Strain = str_replace(Strain, "Strain_", "")) %>% 
  mutate(RR_CI_lab = if_else(Strain==Ref_Strain,
                             "Reference",
                             if_else(is.na(RR_CI_lab),
                                     "n=0",
                                     RR_CI_lab))) %>%
  pivot_wider(names_from = Strain, values_from = RR_CI_lab)

# Recode labels and re-order
est_table_df <- est_table_df %>%
  mutate(Ref_Strain = factor(Ref_Strain, 
                             levels = c("Pre_alpha", "Alpha", "Beta","Delta", "Epsilon", "Eta", "Omicron", "mixed"),
                             labels = c("Pre-alpha", "Alpha", "Beta","Delta", "Epsilon", "Eta", "Omicron", "Mixed"))) %>%
  mutate(Outcome = recode_factor(Outcome, 
                                 # Adverse Birth outcomes
                                 "verylowbirthweight" = "vLBW", 
                                 "lowbirthweight" = "LBW", 
                                 "extremesga" = "eSGA", 
                                 "sga" = "SGA",
                                 "verypreterm" = "vPTB", 
                                 "verypreterm2" = "vPTB (COVID-19 onset <34w)",
                                 "preterm" = "PTB",
                                 "preterm2" = "PTB (COVID-19 onset <37w)",
                                 "fetal_composite_outcome" = "Fetal composite outcome",
                                 # fetal & neonatal morbidity & mortality outcomes:
                                 "stillbirth_28" = "Stillbrith", 
                                 "perinatal_death28" = "Perinatal death", 
                                 "earlyneo_death" = "Early neonatal death",  
                                 "neonatal_death" = "Neonatal death",
                                 "nicu" = "NICU",
                                 # Maternal morbidity & mortality
                                 "prdeath" = "Pregnancy-related death", 
                                 "place_abrupt" = "Placental abruption", 
                                 "preterm_labor" = "Preterm labor", 
                                 "preterm_labor2" = "Preterm labor (COVID-19 onset <37w)",
                                 "haemorrhage" = "Heamorrhage", 
                                 "embolicdz" = "Embolic disease",
                                 "preeclampsia" = "Preeclampsia", 
                                 "eclampsia" = "Eclampsia", 
                                 "pre_or_eclampsia" = "Preeclampsia or eclampsia", 
                                 "hpd_any" = "HDP (diagnosed at any time)",
                                 "hpd_postcovid" = "HDP (diagnosed at or after COVID-19)",
                                 "c_section" = "C-section",
                                 "c_intrap" = "Intrapartum C-section",
                                 "mat_composite_outcome" = "Maternal composite outcome",
                                 # Critical care indicators
                                 "covid_hosp" = "Hospitalization",
                                 "ICUadmit" = "ICU admission",
                                 "critcare" = "Critical care",
                                 "ventilation" = "Ventilation",
                                 "pneumonia" = "Pneumonia")) %>%
  ###############!!!!!!!!!!!!!!!!!!!!!!!!!#################################################
# When threshold is set to ≥70%, manually remove Beta and Epsilon from the line below!! #
#########################################################################################
arrange(Outcome, Ref_Strain) %>%
  select(Outcome, `Reference Strain` = Ref_Strain, `Pre-alpha` = Pre_alpha, Alpha, Beta, Delta, Epsilon, Eta, Omicron, Mixed = mixed)

# Write out est_table_df as a CSV
day_string <- format(Sys.time(), "%Y%m%d")
dir_folder_name <- paste0(day_string, "_dom_strain_", as.character(dom_threshold*100))
dir.create(paste0("data_out/pooled_tables/", dir_folder_name), recursive = TRUE)
write.csv(est_table_df,
          paste0("data_out/pooled_tables/",
                 dir_folder_name, "/Variant_study_pooled_table_", day_string, "_dom_strain_", as.character(dom_threshold*100), ".csv"), row.names = FALSE)

# Create heatmaps, one for each reference strain

outcome_vec <- c("covid_hosp", "ICUadmit", "critcare", "ventilation", "pneumonia",
                 "prdeath", "place_abrupt", "preterm_labor", "preterm_labor2", "haemorrhage", "embolicdz",
                 "preeclampsia", "eclampsia", "pre_or_eclampsia", "hpd_any", "hpd_postcovid", "c_section",
                 "c_intrap","mat_composite_outcome",
                 # fetal outcomes: 
                 "stillbirth_28",  "neonatal_death", "earlyneo_death", 
                 "perinatal_death28", "nicu","lowbirthweight", "verylowbirthweight",
                 "verypreterm", "verypreterm2", "preterm", "preterm2", "sga", "extremesga", "fetal_composite_outcome")
#  "stillbirth_site", "sab_site", "tab_site", "perinatal_deathsite"

# Create a list of plots
plot_list <- list()
# Loop through each reference strain (rs)
for (rs in levels(factor(heatmap_data_df$Ref_Strain))) {
  heatmap_data <- heatmap_data_df %>%
    filter(Ref_Strain==rs)
  
  # Conditional Formatting:
  heatmap_data <- heatmap_data %>% 
    mutate(RR_sig = case_when((RR>1.05 & LCI>1 & UCI>1) ~ "RR>1.05, CI-sig",
                              (RR<0.95 & LCI <1 & UCI<1) ~ "RR<0.95, CI-sig",
                              TRUE ~ "NS"))
  
  heatmap_data$RR_sig <- factor(heatmap_data$RR_sig, 
                                levels=c("RR>1.05, CI-sig", "RR<0.95, CI-sig", 
                                         "NS"),
                                ordered = TRUE)
  
  heatmap_data$Outcome <- factor(heatmap_data$Outcome, 
                                 levels=outcome_vec, 
                                 ordered = TRUE)
  
  myColorScale <- brewer.pal(4,"Set2")
  names(myColorScale) <- levels(heatmap_data$RR_sig)
  levels(heatmap_data$RR_sig)
  
  colScale <- c("RR>1.05, CI-sig" = "orangered", 
                "RR<0.95, CI-sig" = "dodgerblue", 
                "NS" = "gray") 
  
  heatmap_data <- heatmap_data %>%
    mutate(Strain = factor(Strain, 
                           levels = c("Pre_alpha", "Alpha", "Beta", "Delta", "Epsilon", "Eta", "Omicron", "mixed"),
                           labels = c("Pre-alpha", "Alpha", "Beta", "Delta", "Epsilon", "Eta", "Omicron", "Mixed"))) %>%
    mutate(Outcome = recode_factor(Outcome, 
                                   # Adverse Birth outcomes
                                   "verylowbirthweight" = "vLBW", 
                                   "lowbirthweight" = "LBW", 
                                   "extremesga" = "eSGA", 
                                   "sga" = "SGA",
                                   "verypreterm" = "vPTB", 
                                   "verypreterm2" = "vPTB (COVID-19 onset <34w)",
                                   "preterm" = "PTB",
                                   "preterm2" = "PTB (COVID-19 onset <37w)",
                                   "fetal_composite_outcome" = "Fetal composite outcome",
                                   # fetal & neonatal morbidity & mortality outcomes:
                                   "stillbirth_28" = "Stillbrith", 
                                   "perinatal_death28" = "Perinatal death", 
                                   "earlyneo_death" = "Early neonatal death",  
                                   "neonatal_death" = "Neonatal death",
                                   "nicu" = "NICU",
                                   # Maternal morbidity & mortality
                                   "prdeath" = "Pregnancy-related death", 
                                   "place_abrupt" = "Placental abruption", 
                                   "preterm_labor" = "Preterm labor", 
                                   "preterm_labor2" = "Preterm labor (COVID-19 onset <37w)",
                                   "haemorrhage" = "Heamorrhage", 
                                   "embolicdz" = "Embolic disease",
                                   "preeclampsia" = "Preeclampsia", 
                                   "eclampsia" = "Eclampsia", 
                                   "pre_or_eclampsia" = "Preeclampsia or eclampsia", 
                                   "hpd_any" = "HDP (diagnosed at any time)",
                                   "hpd_postcovid" = "HDP (diagnosed at or after COVID-19)",
                                   "c_section" = "C-section",
                                   "c_intrap" = "Intrapartum C-section",
                                   "mat_composite_outcome" = "Maternal composite outcome",
                                   # Critical care indicators
                                   "covid_hosp" = "Hospitalization",
                                   "ICUadmit" = "ICU admission",
                                   "critcare" = "Critical care",
                                   "ventilation" = "Ventilation",
                                   "pneumonia" = "Pneumonia")) 
  # "stillbirth_site","perinatal_deathsite", "sab_site", "tab_site",
  
  hmap <- heatmap_data %>%
    filter(n_study>=1) %>%
    drop_na(Outcome) %>% # Removes NA estimates
    ggplot() +
    geom_tile(color="white", alpha=0.8, aes(x = Strain, y=Outcome, fill=RR_sig))+
    scale_fill_manual(values = colScale) +
    # geom_text(aes(label=RR_CI_lab), size=2) + 
    guides(fill=guide_legend((title = "RR & 95%CI"))) +
    labs(title=paste0("Reference strain: ",rs)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(levels(heatmap_data$Strain)) + # xlim levels of heatmap_df$Strain retains all strains. 
    ylim(levels(heatmap_data$Outcome)) + 
    theme(panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  #legend.position="none")
  
  # Save the plot
  day_string <- format(Sys.time(), "%Y%m%d")
  dir_folder_name <- paste0(day_string,"_dom_strain_",as.character(dom_threshold*100))
  dir.create(paste0("plots/heatmap/",dir_folder_name), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(paste0("plots/heatmap/", dir_folder_name, "/variant_study_heatmap_refstrain_", rs, "_n_study>1.pdf"),
         hmap, width = 8, height = 6, units = "in")
  
  variant_plot <- hmap + theme(legend.position="none")
  plot_list <- append(plot_list, list(variant_plot))
}
# Uncomment to display in RStudi'o:
#
# hmap_layout <- plot_layout(ncol = 2)
# plot_list[[1]] + hmap_layout
