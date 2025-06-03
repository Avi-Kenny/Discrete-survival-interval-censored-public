###########################.
##### Data processing #####
###########################.

chk(1, "Data processing: START")
.t_start <- Sys.time()

set.seed(1)

# Set up data frame to track filtering
dat_log <- data.frame("note"=character(), "num_rows"=integer())
log_note <- function(note, num_rows) {
  dat_log[nrow(dat_log)+1,] <<- list(note, num_rows)
}

# Read in data
dat_prc <- readRDS("../Data/pip_combined_hiv_2025-01-01.rds")
dat_prc %<>% dplyr::ungroup()

log_note("# rows, original", nrow(dat_prc))
log_note("# individuals, original", length(unique(dat_prc$IIntId)))

# Take sample from dataset (for model development)
if (cfg$samp_size!=0) {
  iintids <- unique(dat_prc$IIntId)
  iintids_sample <- sample(iintids, size=cfg$samp_size)
  dat_prc %<>% dplyr::filter(IIntId %in% iintids_sample)
}
log_note("# rows, subsample", nrow(dat_prc))

# Rename columns
dat_prc %<>% dplyr::rename(
  "id" = IIntId,
  "year" = Year
)

# Drop unnecessary columns
dat_prc %<>% subset(select=-c(LocationId, EarliestARTInitDate, HIV_update, age_start,
                              age_end, first_hiv_pos, last_hiv_neg))

# Function to convert dates
convert_date <- memoise(function(date) {
  date <- as.Date(date, format="%Y-%m-%d")
  return(lubridate::year(date))
})

# Filter out obs with sex=="Unknown"
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(sex!="Unknown")
log_note("# rows dropped, sex=='Unknown'", rows_pre-nrow(dat_prc))

# Filter out obs with missing `PIPSA`
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(!is.na(PIPSA))
log_note("# rows dropped, missing `PIPSA`", rows_pre-nrow(dat_prc))

# Filter out PIPSA North
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(PIPSA!="N")
log_note("# rows dropped, PIPSA='N'", rows_pre-nrow(dat_prc))

# Filter dataset based on sex
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(sex==cfg$model_sex)
if (nrow(dat_prc)==0) { stop("`model_sex` incorrectly specified.") }
log_note("# rows after filtering by sex", rows_pre-nrow(dat_prc))

# Misc data wrangling
dat_prc %<>% dplyr::mutate(
  died = ifelse(is.na(died), 0, died),
  dob = convert_date(dob),
  age = year - dob,
  ART_update = ifelse(is.na(ART_update), 0, ART_update)
)

# Filter out children with tests before age 12
rows_pre <- nrow(dat_prc)
children_with_tests <- dplyr::filter(
  dat_prc, age<=13 & ResultDate!=""
)$id
dat_prc %<>% dplyr::filter(!(id %in% children_with_tests))
log_note("# rows dropped, children with tests before age 12",
         rows_pre-nrow(dat_prc))

# Remove all data before 13th birthday
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(age>=13)
log_note("# rows dropped, age<13", rows_pre-nrow(dat_prc))

# Remove all data after age cfg$age_end
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(age<cfg$age_end)
log_note("# rows dropped, age>=cfg$age_end", rows_pre-nrow(dat_prc))

# Remove tests with an "S" result
rows_with_s <- which(dat_prc$HIVResult=="S")
dat_prc[rows_with_s, "ResultDate"] <- NA
dat_prc[rows_with_s, "HIVResult"] <- ""

# Create flag for individuals whose status is positive/known at window_start
dat_pos_at_start <- dat_prc %>%
  dplyr::filter(year<=cfg$w_start) %>%
  mutate(pos=ifelse(is.na(HIVResult), 0, In(HIVResult=="P"))) %>%
  dplyr::group_by(id) %>%
  dplyr::summarize(pos_at_start = max(pos)) %>%
  dplyr::filter(pos_at_start==1)
ids_pos_at_start <- dat_pos_at_start$id
rm(dat_pos_at_start)

# Remove all data prior to window_start
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(year>=cfg$w_start)
log_note("# rows dropped, year<cfg$w_start", rows_pre-nrow(dat_prc))

# Set V=1 if status is known at window_start
rows <- which(
  dat_prc$id %in% ids_pos_at_start & dat_prc$year==cfg$w_start
)
for (row in rows) {
  hiv_res <- dat_prc[row,"HIVResult"]
  if (is.na(hiv_res) || hiv_res!="P") {
    dat_prc[row,"HIVResult"] <- "P"
    dat_prc[row,"ResultDate"] <- as.Date(
      paste0(cfg$w_start, "-01-01"),
      format = "%Y-%m-%d"
    )
  }
}

# Remove all data following window_end
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(year<=cfg$w_end)
log_note("# rows dropped, year>cfg$w_end", rows_pre-nrow(dat_prc))

# Add `first_hiv_pos_dt` and `last_hiv_neg_dt`
dat_prc %<>% dplyr::mutate(
  hiv_pos_dts = ifelse(HIVResult=="P", year, NA),
  hiv_pos_dts = ifelse(is.na(hiv_pos_dts), 9999, hiv_pos_dts),
  hiv_neg_dts = ifelse(HIVResult=="N", year, NA),
  hiv_neg_dts = ifelse(is.na(hiv_neg_dts), 0, hiv_neg_dts),
  art_pos_dts = ifelse(ART_update==1, year, NA),
  art_pos_dts = ifelse(is.na(art_pos_dts), 9999, art_pos_dts)
)
dat_prc %<>% dplyr::mutate(
  first_hiv_pos_dt = min(hiv_pos_dts),
  first_hiv_pos_dt = ifelse(first_hiv_pos_dt==9999, NA, first_hiv_pos_dt),
  last_hiv_neg_dt = max(hiv_neg_dts),
  last_hiv_neg_dt = ifelse(last_hiv_neg_dt==0, NA, last_hiv_neg_dt),
  first_art_pos_dt = min(art_pos_dts),
  .by = id
)
dat_prc[["hiv_pos_dts"]] <- NULL
dat_prc[["hiv_neg_dts"]] <- NULL
dat_prc[["art_pos_dts"]] <- NULL

# Create ART status variable
dat_prc %<>% dplyr::mutate(
  ART_status_new = ifelse(year>=first_art_pos_dt, 1, 0),
  .by = id
)

# Filter out records with a negative test after a positive test
rows_pre <- nrow(dat_prc)
dat_prc %<>% dplyr::filter(
  is.na(first_hiv_pos_dt) | is.na(last_hiv_neg_dt) |
    first_hiv_pos_dt>last_hiv_neg_dt
)
log_note("# rows dropped, NEG test after POS test", rows_pre-nrow(dat_prc))

# Sort dataframe
dat_prc %<>% dplyr::arrange(id,year)

# Drop rows with duplicate person-time
rows_pre <- nrow(dat_prc)
dupe_time_rows <- which(
  dat_prc$id==c(NA,dat_prc$id[1:length(dat_prc$id)-1]) &
    dat_prc$year==c(NA,dat_prc$year[1:length(dat_prc$year)-1])
)
if (length(dupe_time_rows)>0) { dat_prc <- dat_prc[-dupe_time_rows,] }
log_note("# rows dropped, duplicate person-time", rows_pre-nrow(dat_prc))
log_note("# rows, final", nrow(dat_prc))
log_note("# individuals, final", length(unique(dat_prc$id)))
log_note("# deaths, final", sum(dat_prc$died))

# Print data log
print(dat_log)

# Rename columns
dat_prc %<>% dplyr::rename(
  "y" = died,
  "z" = ART_status_new,
  "t_end" = year
)

# Create grouped dataset
dat_grp <- dat_prc %>% dplyr::group_by(id) %>%
  dplyr::summarize(
    count = n(),
    T_minus = last_hiv_neg_dt[1],
    T_plus = first_hiv_pos_dt[1],
    s_i = min(t_end),
    t_i = max(t_end)
  )
dat_grp %<>% dplyr::mutate(
  case = case_when(
    is.na(T_minus) & is.na(T_plus) ~ 1,
    !is.na(T_minus) & is.na(T_plus) ~ 2,
    !is.na(T_minus) & !is.na(T_plus) ~ 3,
    is.na(T_minus) & !is.na(T_plus) ~ 4,
    TRUE ~ 999
  )
)

# Renumber IDs
dat_prc %<>% dplyr::mutate(
  id_orig = id,
  id = as.integer(factor(id))
)

# Set data attributes
attr(dat_prc, "n") <- as.integer(max(dat_prc$id))
attr(dat_prc, "s_i") <- dat_grp$s_i
attr(dat_prc, "t_i") <- dat_grp$t_i

# Generate delta column
delta <- rep(NA, nrow(dat_prc))
for (id in c(1:attr(dat_prc, "n"))) {
  rows_i <- which(dat_prc$id==id)
  delta_i <- g_delta(
    case = dat_grp$case[id],
    s_i = dat_grp$s_i[id],
    t_i = dat_grp$t_i[id],
    T_minus = dat_grp$T_minus[id],
    T_plus = dat_grp$T_plus[id]
  )
  if (length(rows_i)!=length(delta_i)) {
    stop(paste0("Error with computation of delta for ID ", id, "."))
  }
  delta[rows_i] <- delta_i
}
dat_prc$delta <- delta

# Create V (testing) and U (positive/known) indicators
dat_prc %<>% dplyr::mutate(
  v = In(!is.na(dat_prc$HIVResult)),
  u = In(!is.na(first_hiv_pos_dt) & t_end>=first_hiv_pos_dt)
)

# Rescale time variables via scale_time function
dat_prc %<>% dplyr::mutate(
  dob = scale_time(dob, st=cfg$w_start),
  t_end = scale_time(t_end, st=cfg$w_start),
  first_hiv_pos_dt = scale_time(first_hiv_pos_dt, st=cfg$w_start),
  last_hiv_neg_dt = scale_time(last_hiv_neg_dt, st=cfg$w_start)
)
attr(dat_prc, "s_i") <- scale_time(attr(dat_prc, "s_i"),
                                   st=cfg$w_start)
attr(dat_prc, "t_i") <- scale_time(attr(dat_prc, "t_i"),
                                   st=cfg$w_start)

# Create scaled age variable (w_1) and geography covariate (w_2)
# Geography covariate: 1="PIPSA North", 0="PIPSA South"
dat_prc %<>% dplyr::mutate(
  w_1 = scale_age(t_end-dob),
  w_2 = 0
)

# Save datasets for validation
if (cfg$save_dat) {
  saveRDS(dat_prc, paste0("../Data/dat_prc_", substr(cfg$model_sex,1,1), "_",
                          format(Sys.time(), "%Y%m%d"), ".rds"))
}

cols_to_drop <- c(
  "DoB", "dob", "age", "ResultDate", "HIVResult", "hiv_result_fill",
  "VisitDate", "ReceivedHIVTestResult", "CurrentlyOnART", "HadPosHIVResult",
  "first_hiv_pos_dt", "last_hiv_neg_dt", "ART_update", "first_art_pos_dt",
  "sex", "id_orig", "PIPSA", "year_begin", "year_end"
)
for (col in cols_to_drop) { dat_prc[[col]] <- NULL }

# Create transformed dataset object
dat_objs <- transform_dataset(
  dat = dat_prc,
  window_start = cfg$w_start,
  window_end = cfg$w_end
)

saveRDS(dat_prc, paste0("../Data/dat_prc_", cfg$model_sex, ".rds"))
saveRDS(dat_objs, paste0("../Data/dat_objs_", cfg$model_sex, ".rds"))

chk(1, "Data processing: END")
print("Runtime (data processing):")
print(Sys.time()-.t_start)
