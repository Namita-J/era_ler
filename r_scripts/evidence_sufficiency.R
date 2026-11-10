# What evidence do we have? Description of the LER data (outputs/LER_all.csv) by intercropping
# system, geography (incl. ERA agro-ecological zone) and time. Counts are of unique papers,
# because one paper can contribute many LER values.
#
# Input (read only): outputs/LER_all.csv; ERA raw extraction (Data.Out) and master codes, used
#   only to fill missing site metadata and to identify water-management practice codes.
# Output: outputs/evidence_sufficiency/{data,lookups,tables,figures}
# Run from the project root (era_ler/).

suppressMessages({
  library(data.table); library(ggplot2); library(sf); library(rnaturalearth)
  library(jsonlite); library(viridis)
})

es_dir <- file.path("outputs", "evidence_sufficiency")
for (d in c("data", "lookups", "tables", "figures")) dir.create(file.path(es_dir, d), recursive = TRUE, showWarnings = FALSE)
tab <- function(x, name) fwrite(x, file.path(es_dir, "tables", paste0(name, ".csv")))
fig <- function(p, name, w, h) {
  ggsave(file.path(es_dir, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(es_dir, "figures", paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
theme_es <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), legend.position = "bottom")
africa_xlim <- c(-20, 55); africa_ylim <- c(-35, 38)

# ---- 0. Data -----------------------------------------------------------------------------
L <- fread(file.path("outputs", "LER_all.csv"), colClasses = list(character = c("M.Year", "Site.ID", "T.Descrip")))
L[, obs_id := .I]
n_input <- nrow(L)

bundle <- file.path("downloaded_data", "era_agronomy_bundle_2025-05-12")
if (!exists("agronomic_metadata")) agronomic_metadata <- fromJSON(list.files(bundle, "^agronomic_.*json$", full.names = TRUE))
if (!exists("mastercodes")) mastercodes <- fromJSON(list.files(bundle, "^era_master_codes.*json$", full.names = TRUE))

# fill missing site metadata (country, coordinates, author, year) from the raw extraction,
# for rows whose metadata could not be copied from a compiled row (raw-only papers)
raw_site <- as.data.table(agronomic_metadata$Data.Out)[, .(
  Country_raw = Country[1], ISO_raw = ISO.3166.1.alpha.3[1],
  Lat_raw = suppressWarnings(as.numeric(Site.LatD[1])), Lon_raw = suppressWarnings(as.numeric(Site.LonD[1])),
  Author_raw = B.Author.Last[1], Date_raw = suppressWarnings(as.numeric(B.Date[1]))
), by = .(Code = B.Code, Site.ID)]
L <- merge(L, raw_site, by = c("Code", "Site.ID"), all.x = TRUE, sort = FALSE)
L[(is.na(Country) | Country == "") & !is.na(Country_raw), `:=`(Country = Country_raw, ISO.3166.1.alpha.3 = ISO_raw)]
L[is.na(Latitude) & !is.na(Lat_raw), `:=`(Latitude = Lat_raw, Longitude = Lon_raw)]
L[is.na(Author) & !is.na(Author_raw), Author := Author_raw]
L[is.na(Date) & !is.na(Date_raw), Date := Date_raw]
L[, c("Country_raw", "ISO_raw", "Lat_raw", "Lon_raw", "Author_raw", "Date_raw") := NULL]
setorder(L, obs_id)

# regions (UN sub-regions from Natural Earth)
world <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
world$admin[world$admin == "United Republic of Tanzania"] <- "Tanzania"
region_lu <- as.data.table(st_drop_geometry(world))[, .(ISO.3166.1.alpha.3 = adm0_a3, Region = subregion)]
L <- merge(L, region_lu, by = "ISO.3166.1.alpha.3", all.x = TRUE, sort = FALSE)
setorder(L, obs_id)

# agro-ecological zone as recorded in ERA
L[, AEZ := fifelse(is.na(AEZ16simple) | AEZ16simple == "", "Not recorded", AEZ16simple)]

# ---- 1. Papers, locations, years --------------------------------------------------------------
L[, paper := Code]
L[, location := fifelse(!is.na(Latitude) & !is.na(Longitude),
                        paste(sprintf("%.2f", Latitude), sprintf("%.2f", Longitude), sep = ", "), NA_character_)]

# experiment years from M.Year (the experiment period recorded by ERA, not the publication year)
L[, years := lapply(regmatches(M.Year, gregexpr("(19|20)[0-9]{2}", M.Year)), function(v) sort(unique(as.integer(v))))]
L[, n_years_obs := lengths(years)]
L[, year_precision := fcase(n_years_obs == 1, "single year", n_years_obs >= 2, "multi-year average", default = "no calendar year")]
L[, year_start := sapply(years, function(v) if (length(v)) min(v) else NA_integer_)]

# ---- 2. Intercropping systems -----------------------------------------------------------------
synonyms <- c("bean" = "common bean", "beans" = "common bean", "soya bean" = "soybean", "amaranthus" = "amaranth",
              "gliricidia sp." = "gliricidia sp", ".pepper" = "pepper", "gum arabic" = "acacia senegal",
              "acacia" = "acacia senegal", "roselle" = "hibiscus", "millet" = "millet (unspecified)",
              "other millet" = "millet (unspecified)", "banana sweet" = "banana", "arabica" = "coffee (arabica)")
split_crops <- function(x) {
  x <- gsub("mucuna.sunhemp", "mucuna/sunhemp", x, fixed = TRUE)
  v <- tolower(trimws(unlist(strsplit(x, "[-/]"))))
  v <- v[nzchar(v) & !is.na(v) & v != "na"]
  v <- ifelse(v %in% names(synonyms), synonyms[v], v)
  sort(unique(v))
}
crop_type_lu <- data.table(crop = c(
  "maize", "sorghum", "pearl millet", "finger millet", "millet (unspecified)", "rice", "wheat", "barley", "rye", "teff",
  "cowpea", "common bean", "soybean", "groundnut", "pigeon pea", "mung bean", "lablab", "fava bean", "lupin",
  "african yam bean", "butter bean", "lima bean", "velvet bean", "mucuna", "sunhemp", "crotalaria ochrolueca",
  "hairy vetch", "common vetch", "purple vetch", "egyptian clover",
  "cassava", "yam", "potato",
  "amaranth", "calabash", "chili", "pepper", "kale", "melon", "watermelon", "okra", "onion", "pumpkin", "tomato",
  "acacia senegal", "coffee (arabica)", "banana", "cardamom", "gliricidia sp", "grevillea robusta",
  "leucaena leucocephala", "senna siamea", "fodder tree",
  "cotton", "sesame", "hibiscus", "unspecified fodder grass", "natural fallow"),
  crop_type = c(rep("cereal", 10), rep("legume", 20), rep("root/tuber", 3), rep("vegetable", 11),
                rep("tree/perennial", 9), "fibre", "oilseed", "other crop", "forage grass", "fallow"))

L[, crops_product := lapply(Product.Simple, split_crops)]
L[, crops_diversity := lapply(fifelse(is.na(Diversity), "", Diversity), split_crops)]
# a single-crop product label is replaced by the Diversity combination when that names 2+ crops
L[, use_diversity := lengths(crops_product) == 1 & lengths(crops_diversity) >= 2]
L[, crops := fifelse(use_diversity, crops_diversity, crops_product)]
L[, combination_source := fifelse(use_diversity, "Diversity (product label named one crop)", "Product")]
L[, combination := sapply(crops, function(v) paste(v, collapse = " + "))]
unknown_crops <- setdiff(unique(unlist(L$crops)), crop_type_lu$crop)
if (length(unknown_crops)) stop("crops without a type in crop_type_lu: ", paste(unknown_crops, collapse = ", "))

type_of <- function(v) crop_type_lu$crop_type[match(v, crop_type_lu$crop)]
system_group_of <- function(v) {
  t <- type_of(v)
  if (any(t == "tree/perennial")) return("Tree-based (agroforestry / perennial)")
  if (length(v) == 1) return("Unidentified (single-crop label)")
  if (length(v) >= 3) return("Multi-crop (3+ crops)")
  s <- paste(sort(t), collapse = " + ")
  switch(s,
    "cereal + legume" = "Cereal + legume",
    "cereal + cereal" = "Cereal + cereal",
    "legume + root/tuber" = "Root/tuber + legume",
    "cereal + root/tuber" = "Root/tuber + cereal",
    "legume + legume" = "Legume + legume",
    if (any(t == "vegetable")) "Vegetable-based" else "Other")
}
L[, system_group := sapply(crops, system_group_of)]
L[, crop_types := sapply(crops, function(v) paste(sort(unique(type_of(v))), collapse = " + "))]

crosswalk <- unique(L[, .(combination, crop_types, system_group)])[order(system_group, combination)]
fwrite(crop_type_lu[crop %in% unique(unlist(L$crops))][order(crop_type, crop)], file.path(es_dir, "lookups", "crop_type_lookup.csv"))
fwrite(crosswalk, file.path(es_dir, "lookups", "system_group_crosswalk.csv"))
group_levels <- c("Cereal + legume", "Cereal + cereal", "Root/tuber + legume", "Root/tuber + cereal", "Legume + legume",
                  "Multi-crop (3+ crops)", "Vegetable-based", "Tree-based (agroforestry / perennial)", "Other",
                  "Unidentified (single-crop label)")
L[, system_group := factor(system_group, levels = intersect(group_levels, unique(system_group)))]

# ---- 3. Water-management / irrigation information ------------------------------------------
# Treatment fields describe the LER observation itself for rows from the compiled table and rows
# calculated from yields; for reported_raw rows they were copied from another row of the paper,
# so they are not used for those rows.
water_codes <- as.data.table(mastercodes$prac)[Theme %in% c("Blue Water Management", "Green Water Management"), Code]
L[, obs_specific_treatment := LER.Source != "reported_raw"]
tcols <- paste0("T", 1:13)
L[, t_water := apply(as.matrix(.SD), 1, function(r) paste(sort(intersect(r, water_codes)), collapse = " ")), .SDcols = tcols]
L[, plist_water := sapply(strsplit(fifelse(is.na(plist), "", plist), "-"), function(v) any(v %in% water_codes))]

# ---- 4. Per-paper summary ----------------------------------------------------------------------
per_paper <- L[, .(
  Author = Author[1], Year_published = Date[1], Country = paste(sort(unique(na.omit(Country))), collapse = "; "),
  n_obs = .N, n_sites = uniqueN(Site.ID),
  n_calendar_years = uniqueN(unlist(years)),
  n_single_years = uniqueN(year_start[year_precision == "single year"]),   # distinct years reported separately
  system_groups = paste(sort(unique(as.character(system_group))), collapse = "; ")
), by = paper]

# ---- 5. Systems ----------------------------------------------------------------------------
group_tab <- L[, .(combinations = uniqueN(combination), papers = uniqueN(paper), countries = uniqueN(na.omit(Country)),
                   LER_values = .N), by = system_group][order(system_group)]
combo_tab <- L[, .(papers = uniqueN(paper), LER_values = .N), by = .(combination, system_group)][order(-papers, -LER_values)]
tab(group_tab, "02_system_groups"); tab(combo_tab, "02_crop_combinations")

# ---- 6. Geography -------------------------------------------------------------------------
region_tab <- L[!is.na(Region), .(countries = uniqueN(Country), papers = uniqueN(paper), LER_values = .N), by = Region][order(-papers)]
country_tab <- L[!is.na(Country) & Country != "", .(papers = uniqueN(paper), LER_values = .N), by = .(Region, Country)][order(-papers)]
aez_tab <- L[, .(papers = uniqueN(paper), countries = uniqueN(na.omit(Country)), LER_values = .N), by = AEZ][order(AEZ == "Not recorded", -papers)]
tab(region_tab, "03_papers_by_region"); tab(country_tab, "03_papers_by_country"); tab(aez_tab, "03_papers_by_aez")

base_map <- function() ggplot() + geom_sf(data = world, fill = "gray95", color = "white", linewidth = 0.2) +
  coord_sf(xlim = africa_xlim, ylim = africa_ylim, expand = FALSE) + theme_es +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())

loc_tab <- L[!is.na(location), .(papers = uniqueN(paper), Latitude = Latitude[1], Longitude = Longitude[1]), by = location]
p_loc <- base_map() +
  geom_point(data = loc_tab, aes(Longitude, Latitude, size = papers), shape = 21, fill = "#2a6f97", color = "black", alpha = 0.7) +
  scale_size_area(max_size = 6, name = "Papers") +
  labs(title = "Study locations", subtitle = paste0(nrow(loc_tab), " locations"))
fig(p_loc, "03a_study_locations", 8, 8)

country_papers <- merge(world, L[!is.na(ISO.3166.1.alpha.3), .(papers = uniqueN(paper)), by = .(adm0_a3 = ISO.3166.1.alpha.3)], by = "adm0_a3", all.x = TRUE)
p_country <- ggplot() + geom_sf(data = country_papers, aes(fill = papers), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "mako", direction = -1, na.value = "gray95", name = "Papers") +
  coord_sf(xlim = africa_xlim, ylim = africa_ylim, expand = FALSE) + theme_es +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
  labs(title = "Papers per country")
fig(p_country, "03b_country_papers", 8, 8)

L[!is.na(location), `:=`(cell_lat = floor(Latitude) + 0.5, cell_lon = floor(Longitude) + 0.5)]
grid_tab <- L[!is.na(location), .(papers = uniqueN(paper)), by = .(cell_lat, cell_lon)]
p_grid <- base_map() +
  geom_tile(data = grid_tab, aes(cell_lon, cell_lat, fill = papers), width = 1, height = 1, color = "gray30", linewidth = 0.1) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Papers") + labs(title = "Papers per 1-degree cell")
fig(p_grid, "03c_grid_papers", 8, 8)

aez_order <- aez_tab$AEZ
p_aez <- ggplot(aez_tab, aes(papers, factor(AEZ, levels = rev(aez_order)))) + geom_col(fill = "#2a6f97", width = 0.7) +
  geom_text(aes(label = papers), hjust = -0.3, size = 3.5) + scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Papers", y = NULL, title = "Papers per agro-ecological zone") + theme_es
fig(p_aez, "03d_aez_papers", 8, 8)

# ---- 7. System x geography -----------------------------------------------------------------
sg_region <- L[!is.na(Region), .(papers = uniqueN(paper)), by = .(system_group, Region)]
sg_aez <- L[, .(papers = uniqueN(paper)), by = .(system_group, AEZ)]
sg_country <- L[!is.na(Country) & Country != "", .(papers = uniqueN(paper)), by = .(system_group, Country)]
tab(sg_region, "04_system_by_region"); tab(sg_aez, "04_system_by_aez")
heat <- function(d, xvar, title, xlev = NULL) {
  d <- copy(d); d[, x := get(xvar)]
  if (is.null(xlev)) xlev <- d[, .(t = sum(papers)), by = x][order(-t), x]
  ggplot(d, aes(factor(x, levels = xlev), system_group, fill = papers)) + geom_tile(color = "white") +
    geom_text(aes(label = papers, color = papers > 0.6 * max(papers)), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.1, end = 0.85, name = "Papers") +
    scale_y_discrete(limits = rev) + labs(x = NULL, y = NULL, title = title) + theme_es +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}
fig(heat(sg_region, "Region", "System group x region: papers"), "04a_system_region_papers", 9, 6)
fig(heat(sg_aez, "AEZ", "System group x agro-ecological zone: papers", aez_order), "04b_system_aez_papers", 9, 6)

# geographic concentration of each system group (share of paper-country pairs)
conc <- sg_country[, .(papers = sum(papers), top_country = Country[which.max(papers)],
                       top_country_share_pct = round(100 * max(papers) / sum(papers))), by = system_group]

# ---- 8. Time -------------------------------------------------------------------------------
year_tab <- L[n_years_obs > 0, .(year = unlist(years)), by = .(obs_id, paper)][, .(papers = uniqueN(paper)), by = year][order(year)]
tab(year_tab, "05_papers_by_year")
p_year <- ggplot(year_tab, aes(year, papers)) + geom_col(fill = "#2a6f97") +
  labs(x = "Experiment year", y = "Papers", title = "Papers by experiment year",
       subtitle = "A paper covering several years is counted in each year") + theme_es
fig(p_year, "05a_papers_by_year", 10, 5)

# number of calendar years each paper reports LER for (years averaged together count as covered)
years_per_paper <- per_paper[, .(papers = .N), by = .(n_years = fifelse(n_calendar_years == 0, "No year", as.character(n_calendar_years)))]
yr_lev <- c("No year", as.character(sort(unique(per_paper[n_calendar_years > 0, n_calendar_years]))))
years_per_paper[, n_years := factor(n_years, levels = yr_lev)]
setorder(years_per_paper, n_years)
tab(years_per_paper, "05_years_per_paper")
p_ypp <- ggplot(years_per_paper, aes(n_years, papers)) + geom_col(fill = "#2a6f97", width = 0.7) +
  geom_text(aes(label = papers), vjust = -0.4, size = 3.5) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of years reported", y = "Papers", title = "Number of years each paper reports LER for") + theme_es
fig(p_ypp, "05b_years_per_paper", 8, 5)

# ---- 9. Papers with repeated conditions ---------------------------------------------------------
season_multi <- L[year_precision == "single year", .(n_seasons = uniqueN(M.Year)), by = .(paper, year = year_start)][, .(multi_season = any(n_seasons > 1)), by = paper]
water <- L[, .(
  irrigated_and_rainfed_obs = any(obs_specific_treatment & Irrigation.T %in% TRUE) & any(obs_specific_treatment & Irrigation.T %in% FALSE),
  irrigation_contrast_in_comparison = any(obs_specific_treatment & !is.na(Irrigation.C) & !is.na(Irrigation.T) & Irrigation.C != Irrigation.T),
  water_mgmt_treatment = any(obs_specific_treatment & plist_water) | uniqueN(t_water[obs_specific_treatment]) > 1
), by = paper]
cand <- merge(per_paper, season_multi, by = "paper", all.x = TRUE)
cand <- merge(cand, water, by = "paper")
cand[is.na(multi_season), multi_season := FALSE]
cand[, `:=`(multi_year = n_single_years >= 2, multi_site = n_sites >= 2)]
cand[, water_contrast := irrigated_and_rainfed_obs | irrigation_contrast_in_comparison | water_mgmt_treatment]
cand[, repeated_conditions := multi_year | multi_season | water_contrast]
setorder(cand, -repeated_conditions, -n_single_years)
tab(cand, "06_repeated_conditions_by_paper")
cand_summary <- data.table(
  feature = c("2+ years reported separately", "2+ seasons within a year", "2+ sites",
              "Water / irrigation contrast", "Any repeated condition (years, seasons or water)"),
  papers = c(sum(cand$multi_year), sum(cand$multi_season), sum(cand$multi_site), sum(cand$water_contrast), sum(cand$repeated_conditions)))
cand_summary[, pct_of_papers := round(100 * papers / nrow(cand), 1)]
tab(cand_summary, "06_repeated_conditions_summary")

# ---- 10. Derived dataset ----------------------------------------------------------------------
derived <- L[, .(obs_id, paper, Author, Date, Country, ISO.3166.1.alpha.3, Region, AEZ, Site.ID, Latitude, Longitude,
                 location, M.Year, year_precision, years = sapply(years, paste, collapse = ";"), year_start,
                 Product.Simple, Diversity, combination, combination_source, crop_types, system_group, T.Descrip, LER = MeanT, LER.Source,
                 Irrigation.C, Irrigation.T, water_practices_T = t_water, water_practice_in_plist = plist_water)]
fwrite(derived, file.path(es_dir, "data", "LER_evidence_structure.csv"))
stopifnot(nrow(derived) == n_input)
