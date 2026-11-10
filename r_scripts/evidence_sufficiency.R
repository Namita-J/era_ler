# Evidence-sufficiency assessment of the clean LER dataset
#
# Structural, geographic, temporal and intercropping-system assessment of the independent
# evidence in outputs/LER_all.csv. NO environmental covariates (SPEI, rainfall, temperature,
# soils, AEZ) are used here.
#
# Input (read only): outputs/LER_all.csv; ERA raw extraction (Data.Out) and master codes, used
#   only to fill missing site metadata and to identify water-management practice codes.
# Output: outputs/evidence_sufficiency/{data,lookups,tables,figures}
#
# Units of evidence (kept distinct throughout):
#   paper            = ERA study code (Code)
#   paper-site       = Code x Site.ID (sites as recorded in ERA; one location can hold several,
#                      e.g. fertility levels recorded as separate sites)
#   location         = unique geographic coordinates (Latitude/Longitude rounded to 0.01 deg),
#                      pooled across papers
#   study year       = calendar year(s) parsed from M.Year (the experiment year, not publication year)
#   site-year        = Code x Site.ID x M.Year as recorded (a multi-year average is ONE site-year)
#   system           = order-independent crop combination; system group = functional grouping
#   observation      = one LER row
#
# Assumptions are listed in outputs/evidence_sufficiency/tables/assumptions.csv.
# Run from the project root (era_ler/).

suppressMessages({
  library(data.table); library(ggplot2); library(sf); library(rnaturalearth)
  library(patchwork); library(jsonlite); library(viridis)
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
fill_log <- L[(is.na(Country) | Country == "") | is.na(Latitude) | is.na(Author) | is.na(Date),
              .(rows = .N, country_filled = sum((is.na(Country) | Country == "") & !is.na(Country_raw)),
                coords_filled = sum(is.na(Latitude) & !is.na(Lat_raw)),
                author_filled = sum(is.na(Author) & !is.na(Author_raw))), by = .(Code, Site.ID, LER.Source)]
L[(is.na(Country) | Country == "") & !is.na(Country_raw), `:=`(Country = Country_raw, ISO.3166.1.alpha.3 = ISO_raw)]
L[is.na(Latitude) & !is.na(Lat_raw), `:=`(Latitude = Lat_raw, Longitude = Lon_raw)]
L[is.na(Author) & !is.na(Author_raw), Author := Author_raw]
L[is.na(Date) & !is.na(Date_raw), Date := Date_raw]
L[, c("Country_raw", "ISO_raw", "Lat_raw", "Lon_raw", "Author_raw", "Date_raw") := NULL]
setorder(L, obs_id)
tab(fill_log, "metadata_filled_from_raw")

# regions (UN sub-regions from Natural Earth)
world <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
world$admin[world$admin == "United Republic of Tanzania"] <- "Tanzania"
region_lu <- as.data.table(st_drop_geometry(world))[, .(ISO.3166.1.alpha.3 = adm0_a3, Region = subregion)]
L <- merge(L, region_lu, by = "ISO.3166.1.alpha.3", all.x = TRUE, sort = FALSE)
setorder(L, obs_id)

# ---- 1. Units of evidence -------------------------------------------------------------------
L[, paper := Code]
L[, paper_site := paste(Code, Site.ID, sep = " | ")]
L[, site_aggregated := grepl("..", Site.ID, fixed = TRUE)]          # several sites pooled in one record
L[, location := fifelse(!is.na(Latitude) & !is.na(Longitude),
                        paste(sprintf("%.2f", Latitude), sprintf("%.2f", Longitude), sep = ", "), NA_character_)]
L[, site_year := paste(Code, Site.ID, M.Year, sep = " | ")]

# experiment years from M.Year (the experiment period recorded by ERA, not the publication year)
L[, years := lapply(regmatches(M.Year, gregexpr("(19|20)[0-9]{2}", M.Year)), function(v) sort(unique(as.integer(v))))]
L[, n_years_obs := lengths(years)]
L[, year_precision := fcase(n_years_obs == 1, "single year", n_years_obs >= 2, "multi-year average", default = "no calendar year")]
L[, year_start := sapply(years, function(v) if (length(v)) min(v) else NA_integer_)]
L[, decade := fifelse(is.na(year_start), NA_character_, paste0(floor(year_start / 10) * 10, "s"))]

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
# Treatment fields (practice codes, irrigation) describe the LER observation itself for rows
# from the compiled table and for rows calculated from yields; for reported_raw rows they were
# copied from another row of the same paper, so only the treatment text is used for those.
water_codes <- as.data.table(mastercodes$prac)[Theme %in% c("Blue Water Management", "Green Water Management"), Code]
L[, obs_specific_treatment := LER.Source != "reported_raw"]
tcols <- paste0("T", 1:13)
L[, t_water := apply(as.matrix(.SD), 1, function(r) paste(sort(intersect(r, water_codes)), collapse = " ")), .SDcols = tcols]
L[, plist_water := sapply(strsplit(fifelse(is.na(plist), "", plist), "-"), function(v) any(v %in% water_codes))]
water_text <- "irrigat|rainfed|rain-fed|rain fed|deficit|supplement|water harvest|zai|tied ridg|drip|full water|water regime|water stress"
L[, water_text_cue := grepl(water_text, T.Descrip, ignore.case = TRUE)]

# ---- 4. Evidence-structure summary -----------------------------------------------------------
per_paper <- L[, .(
  Author = Author[1], Year_published = Date[1], Country = paste(sort(unique(na.omit(Country))), collapse = "; "),
  Region = paste(sort(unique(na.omit(Region))), collapse = "; "),
  LER_source = paste(sort(unique(LER.Source)), collapse = "; "),
  n_obs = .N, n_paper_sites = uniqueN(Site.ID), n_locations = uniqueN(na.omit(location)),
  n_site_years = uniqueN(site_year),
  n_calendar_years = uniqueN(unlist(years)),
  n_single_years = uniqueN(year_start[year_precision == "single year"]),   # distinct calendar years observed separately
  any_multi_year_average = any(year_precision == "multi-year average"),
  any_no_calendar_year = any(year_precision == "no calendar year"),
  first_year = suppressWarnings(as.numeric(min(unlist(years)))), last_year = suppressWarnings(as.numeric(max(unlist(years)))),
  n_systems = uniqueN(combination), n_system_groups = uniqueN(system_group),
  system_groups = paste(sort(unique(as.character(system_group))), collapse = "; ")
), by = paper]
per_paper[!is.finite(first_year), `:=`(first_year = NA, last_year = NA)]
setorder(per_paper, -n_obs)

top_share <- function(x, k) round(100 * sum(sort(x, decreasing = TRUE)[seq_len(min(k, length(x)))]) / sum(x), 1)
q <- quantile(per_paper$n_obs, c(.25, .5, .75))
summary_tab <- data.table(
  metric = c("Unique papers", "LER observations", "Paper-sites (paper x Site.ID)", "Geographic locations (unique coordinates, 0.01 deg)",
             "Paper-sites without coordinates", "Paper-sites recording several pooled sites", "Countries",
             "Unique study (calendar) years", "Site-years (paper x site x M.Year)",
             "Site-years that are multi-year averages", "Site-years with no calendar year",
             "LER observations per paper: median", "LER observations per paper: IQR", "LER observations per paper: range",
             "Papers with one paper-site / several", "Papers with one location / several",
             "Papers covering one calendar year / several", "Papers with one system (combination) / several",
             "Papers with one system group / several",
             "% of observations from the top 5 papers", "% of observations from the top 10 papers",
             "% of site-years from the top 5 papers (ranked by site-years)", "% of site-years from the top 10 papers (ranked by site-years)"),
  value = c(uniqueN(L$paper), nrow(L), uniqueN(L$paper_site), uniqueN(na.omit(L$location)),
            uniqueN(L[is.na(location), paper_site]), uniqueN(L[site_aggregated == TRUE, paper_site]), uniqueN(na.omit(L$Country)),
            uniqueN(unlist(L$years)), uniqueN(L$site_year),
            uniqueN(L[year_precision == "multi-year average", site_year]), uniqueN(L[year_precision == "no calendar year", site_year]),
            q[[2]], paste0(q[[1]], " - ", q[[3]]), paste0(min(per_paper$n_obs), " - ", max(per_paper$n_obs)),
            paste0(sum(per_paper$n_paper_sites == 1), " / ", sum(per_paper$n_paper_sites > 1)),
            paste0(sum(per_paper$n_locations <= 1), " / ", sum(per_paper$n_locations > 1)),
            paste0(sum(per_paper$n_calendar_years <= 1), " / ", sum(per_paper$n_calendar_years > 1)),
            paste0(sum(per_paper$n_systems == 1), " / ", sum(per_paper$n_systems > 1)),
            paste0(sum(per_paper$n_system_groups == 1), " / ", sum(per_paper$n_system_groups > 1)),
            top_share(per_paper$n_obs, 5), top_share(per_paper$n_obs, 10),
            top_share(per_paper$n_site_years, 5), top_share(per_paper$n_site_years, 10)))
summary_tab[, value := as.character(value)]
tab(summary_tab, "01_evidence_structure_summary")
tab(per_paper, "01_per_paper_evidence")

# missingness / exclusions (nothing is dropped; rows lacking a field are excluded only from
# the analyses that need that field)
missing_tab <- data.table(
  field = c("Country", "Coordinates", "Calendar experiment year", "System (crop combination identified)"),
  rows_missing = c(L[is.na(Country) | Country == "", .N], L[is.na(location), .N], L[year_precision == "no calendar year", .N],
                   L[system_group == "Unidentified (single-crop label)", .N]),
  papers_affected = c(L[is.na(Country) | Country == "", uniqueN(paper)], L[is.na(location), uniqueN(paper)],
                      L[year_precision == "no calendar year", uniqueN(paper)], L[system_group == "Unidentified (single-crop label)", uniqueN(paper)]),
  excluded_from = c("country / region tables and maps", "location map and grid maps", "year / decade figures",
                    "kept as its own group in system analyses"))
tab(missing_tab, "00_missingness_and_exclusions")

# ---- 5. Systems ----------------------------------------------------------------------------
combo_tab <- L[, .(papers = uniqueN(paper), site_years = uniqueN(site_year), observations = .N,
                   combination_from_diversity = sum(use_diversity)), by = .(combination, system_group)][order(-papers, -observations)]
group_tab <- L[, .(combinations = uniqueN(combination), papers = uniqueN(paper), site_years = uniqueN(site_year),
                   locations = uniqueN(na.omit(location)), countries = uniqueN(na.omit(Country)), observations = .N), by = system_group][order(system_group)]
tab(combo_tab, "02_crop_combinations")
tab(group_tab, "02_system_groups")

# ---- 6. Geography -------------------------------------------------------------------------
country_tab <- L[!is.na(Country) & Country != "", .(papers = uniqueN(paper), paper_sites = uniqueN(paper_site),
                   locations = uniqueN(na.omit(location)), site_years = uniqueN(site_year), observations = .N,
                   system_groups = uniqueN(system_group)), by = .(Region, Country)][order(-papers)]
region_tab <- L[!is.na(Region), .(countries = uniqueN(Country), papers = uniqueN(paper), paper_sites = uniqueN(paper_site),
                  locations = uniqueN(na.omit(location)), site_years = uniqueN(site_year), observations = .N,
                  system_groups = uniqueN(system_group)), by = Region][order(-papers)]
tab(country_tab, "03_country_evidence"); tab(region_tab, "03_region_evidence")

base_map <- function() ggplot() + geom_sf(data = world, fill = "gray95", color = "white", linewidth = 0.2) +
  coord_sf(xlim = africa_xlim, ylim = africa_ylim, expand = FALSE) + theme_es +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())

loc_tab <- L[!is.na(location), .(papers = uniqueN(paper), site_years = uniqueN(site_year), Latitude = Latitude[1], Longitude = Longitude[1]), by = location]
p_loc <- base_map() +
  geom_point(data = loc_tab, aes(Longitude, Latitude, size = site_years, fill = papers), shape = 21, color = "black", alpha = 0.8) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Papers at location") +
  scale_size_area(max_size = 6, name = "Site-years") +
  labs(title = "A. Unique study locations", subtitle = paste0(nrow(loc_tab), " locations; ", L[is.na(location), uniqueN(paper_site)], " paper-sites lack coordinates"))
fig(p_loc, "03a_study_locations", 8, 8)

country_papers <- merge(world, L[!is.na(ISO.3166.1.alpha.3), .(papers = uniqueN(paper)), by = .(adm0_a3 = ISO.3166.1.alpha.3)], by = "adm0_a3", all.x = TRUE)
p_country <- ggplot() + geom_sf(data = country_papers, aes(fill = papers), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "mako", direction = -1, na.value = "gray95", name = "Unique papers") +
  coord_sf(xlim = africa_xlim, ylim = africa_ylim, expand = FALSE) + theme_es +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) +
  labs(title = "B. Unique papers per country")
fig(p_country, "03b_country_papers", 8, 8)

# 1 x 1 degree grid (coordinates are mostly precise: see coordinate_precision.csv)
tab(L[!is.na(location), .(rows = .N, papers = uniqueN(paper)), by = .(buffer_m = cut(Buffer, c(0, 1000, 10000, 50000, Inf), dig.lab = 6))], "03_coordinate_precision")
L[!is.na(location), `:=`(cell_lat = floor(Latitude) + 0.5, cell_lon = floor(Longitude) + 0.5)]
grid_tab <- L[!is.na(location), .(papers = uniqueN(paper), site_years = uniqueN(site_year), observations = .N), by = .(cell_lat, cell_lon)]
tab(grid_tab[order(-papers)], "03_grid_1deg_evidence")
grid_map <- function(var, title, lab) base_map() +
  geom_tile(data = grid_tab, aes(cell_lon, cell_lat, fill = .data[[var]]), width = 1, height = 1, color = "gray30", linewidth = 0.1) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = lab) + labs(title = title)
fig(grid_map("papers", "C. Unique papers per 1-degree cell", "Papers"), "03c_grid_papers", 8, 8)
fig(grid_map("site_years", "D. Site-years per 1-degree cell", "Site-years"), "03d_grid_site_years", 8, 8)

# ---- 7. System x geography -----------------------------------------------------------------
sg_region <- L[!is.na(Region), .(papers = uniqueN(paper), site_years = uniqueN(site_year)), by = .(system_group, Region)]
sg_country <- L[!is.na(Country) & Country != "", .(papers = uniqueN(paper), site_years = uniqueN(site_year)), by = .(system_group, Country)]
tab(sg_region, "04_system_by_region"); tab(sg_country, "04_system_by_country")
heat <- function(d, xvar, val, title) {
  d <- copy(d); d[, x := get(xvar)]
  ord <- d[, .(t = sum(get(val))), by = x][order(-t), x]
  ggplot(d, aes(factor(x, levels = ord), system_group, fill = .data[[val]])) + geom_tile(color = "white") +
    geom_text(aes(label = .data[[val]], color = .data[[val]] > 0.6 * max(.data[[val]])), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.1, end = 0.85, name = gsub("_", " ", val)) +
    scale_y_discrete(limits = rev) + labs(x = NULL, y = NULL, title = title) + theme_es +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}
fig(heat(sg_region, "Region", "papers", "System group x region: unique papers"), "04a_system_region_papers", 9, 6)
fig(heat(sg_region, "Region", "site_years", "System group x region: site-years"), "04b_system_region_site_years", 9, 6)
fig(heat(sg_country, "Country", "papers", "System group x country: unique papers"), "04c_system_country_papers", 13, 6.5)

# geographic concentration of each system group
conc <- sg_country[, .(papers = sum(papers), top_country = Country[which.max(papers)],
                       top_country_share_pct = round(100 * max(papers) / sum(papers))), by = system_group]
conc <- merge(conc, sg_region[, .(top_region = Region[which.max(papers)], top_region_share_pct = round(100 * max(papers) / sum(papers))), by = system_group], by = "system_group")
# note: a paper spanning two countries counts in both, so shares are of paper-country pairs
tab(conc[order(-papers)], "04_system_geographic_concentration")

# ---- 8. Time -------------------------------------------------------------------------------
yr_long <- L[n_years_obs > 0, .(year = unlist(years)), by = .(obs_id, paper, site_year, system_group, year_precision)]
year_tab <- yr_long[, .(papers = uniqueN(paper), site_years = uniqueN(site_year),
                        site_years_single_year = uniqueN(site_year[year_precision == "single year"])), by = year][order(year)]
decade_tab <- L[!is.na(decade), .(papers = uniqueN(paper), site_years = uniqueN(site_year), observations = .N), by = decade][order(decade)]
sg_decade <- L[!is.na(decade), .(papers = uniqueN(paper), site_years = uniqueN(site_year)), by = .(system_group, decade)]
tab(year_tab, "05_evidence_by_year"); tab(decade_tab, "05_evidence_by_decade"); tab(sg_decade, "05_system_by_decade")
precision_tab <- per_paper[, .(papers = .N), by = .(
  year_information = fcase(n_single_years > 0 & !any_multi_year_average & !any_no_calendar_year, "single-year observations only",
                           n_single_years > 0, "some single-year observations",
                           any_multi_year_average & !any_no_calendar_year, "multi-year averages only",
                           default = "no calendar year")
)][, pct := round(100 * papers / sum(papers), 1)][]
years_per_study <- per_paper[, .(papers = .N), by = .(calendar_years_covered = fcase(n_calendar_years == 0, "0 (no calendar year)",
  n_calendar_years == 1, "1", n_calendar_years == 2, "2", default = "3 or more"))][order(calendar_years_covered)]
years_separately_observed <- per_paper[, .(papers = .N), by = .(separately_observed_single_years = fcase(n_single_years == 0, "0",
  n_single_years == 1, "1", n_single_years == 2, "2", default = "3 or more"))][order(separately_observed_single_years)]
tab(precision_tab, "05_year_precision_by_paper"); tab(years_per_study, "05_years_per_study"); tab(years_separately_observed, "05_separately_observed_years_per_study")

p_year <- ggplot(melt(year_tab[, .(year, papers, site_years)], id.vars = "year"), aes(year, value, fill = variable)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c(papers = "#2a6f97", site_years = "#a9d6e5"), labels = c("Unique papers", "Site-years"), name = NULL) +
  labs(x = "Experiment year", y = "Count", title = "Evidence by experiment year",
       subtitle = "A paper or site-year spanning several years is counted in each year it covers") + theme_es
fig(p_year, "05a_evidence_by_year", 10, 5)
fig(heat(sg_decade, "decade", "papers", "System group x decade (experiment start): unique papers") +
      scale_x_discrete(limits = sort(unique(sg_decade$decade))), "05b_system_decade_papers", 8, 6)

# ---- 9. Repeated-condition / adaptation-analysis candidates ---------------------------------
season_multi <- L[year_precision == "single year", .(n_seasons = uniqueN(M.Year)), by = .(paper, year = year_start)][, .(multi_season = any(n_seasons > 1)), by = paper]
water <- L[, .(
  irrigated_and_rainfed_obs = any(obs_specific_treatment & Irrigation.T %in% TRUE) & any(obs_specific_treatment & Irrigation.T %in% FALSE),
  irrigation_contrast_in_comparison = any(obs_specific_treatment & !is.na(Irrigation.C) & !is.na(Irrigation.T) & Irrigation.C != Irrigation.T),
  water_mgmt_treatment = any(obs_specific_treatment & plist_water) | uniqueN(t_water[obs_specific_treatment]) > 1,
  water_text_cue_only = any(water_text_cue)
), by = paper]
cand <- merge(per_paper[, .(paper, Author, Year_published, Country, LER_source, n_obs, n_paper_sites, n_locations, n_site_years,
                            n_calendar_years, n_single_years, system_groups)], season_multi, by = "paper", all.x = TRUE)
cand <- merge(cand, water, by = "paper")
cand[is.na(multi_season), multi_season := FALSE]
cand[, `:=`(multi_year = n_single_years >= 2, multi_year_average_only = n_calendar_years >= 2 & n_single_years < 2,
            multi_site = n_paper_sites >= 2, multi_location = n_locations >= 2)]
cand[, water_contrast := irrigated_and_rainfed_obs | irrigation_contrast_in_comparison | water_mgmt_treatment]
cand[, water_text_cue_only := water_text_cue_only & !water_contrast]
# repeated-condition comparison within a study: the same system observed in separately recorded
# years or seasons, or under a recorded water/irrigation contrast
cand[, within_study_repeated := multi_year | multi_season | water_contrast]
setcolorder(cand, c("paper", "Author", "Year_published", "Country", "n_obs", "n_site_years", "multi_year", "multi_year_average_only",
                    "multi_site", "multi_location", "multi_season", "irrigated_and_rainfed_obs", "irrigation_contrast_in_comparison",
                    "water_mgmt_treatment", "water_contrast", "water_text_cue_only", "within_study_repeated"))
setorder(cand, -within_study_repeated, -n_site_years)
tab(cand, "06_adaptation_candidates_by_paper")
cand_summary <- data.table(
  feature = c("Multi-year (2+ separately observed years)", "Multi-year averages only (years pooled)", "Multi-site (2+ paper-sites)",
              "Multi-location (2+ distinct coordinates)", "Multi-season (2+ seasons within a year)",
              "Irrigated and rainfed observations", "Irrigation contrast within a comparison", "Water-management treatment",
              "Any water / irrigation contrast", "Water wording in treatment text only (unverified)",
              "Potential within-study repeated-condition comparison"),
  papers = c(sum(cand$multi_year), sum(cand$multi_year_average_only), sum(cand$multi_site), sum(cand$multi_location), sum(cand$multi_season),
             sum(cand$irrigated_and_rainfed_obs), sum(cand$irrigation_contrast_in_comparison), sum(cand$water_mgmt_treatment),
             sum(cand$water_contrast), sum(cand$water_text_cue_only), sum(cand$within_study_repeated)))
cand_summary[, pct_of_papers := round(100 * papers / nrow(cand), 1)]
tab(cand_summary, "06_adaptation_candidates_summary")

# ---- 10. Pseudo-replication -----------------------------------------------------------------
pp <- copy(per_paper)[, paper := factor(paper, levels = paper)]
p_obs <- ggplot(pp, aes(paper, n_obs)) + geom_col(fill = "#2a6f97", width = 0.85) +
  labs(x = "Paper (ordered by LER observations)", y = "LER observations",
       title = sprintf("%d LER observations come from %d papers", nrow(L), nrow(per_paper)),
       subtitle = sprintf("Top 5 papers = %s%% of observations; top 10 = %s%%", top_share(per_paper$n_obs, 5), top_share(per_paper$n_obs, 10))) +
  theme_es + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
p_sy <- ggplot(pp, aes(paper, n_site_years)) + geom_col(fill = "#a9d6e5", color = "#2a6f97", linewidth = 0.1, width = 0.85) +
  labs(x = "Paper (same order)", y = "Site-years", title = sprintf("...and %d site-years", uniqueN(L$site_year)),
       subtitle = sprintf("Median %s site-years per paper", median(per_paper$n_site_years))) +
  theme_es + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
fig(p_obs / p_sy, "07_pseudo_replication", 11, 8)

# ---- 11. Four-panel evidence-sufficiency figure -------------------------------------------------
gl <- melt(group_tab[, .(system_group, Papers = papers, `Site-years` = site_years)], id.vars = "system_group")
p_b <- ggplot(gl, aes(value, system_group, fill = variable)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = value), position = position_dodge(width = 0.8), hjust = -0.15, size = 2.8) +
  scale_fill_manual(values = c(Papers = "#2a6f97", `Site-years` = "#a9d6e5"), name = NULL) +
  scale_y_discrete(limits = rev) + scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Count", y = NULL, title = "B. Independent evidence by system group") + theme_es
p_d <- ggplot(pp, aes(paper, n_obs)) + geom_col(fill = "#2a6f97", width = 0.85) +
  labs(x = "Papers (ordered)", y = "LER observations", title = "D. LER observations per paper",
       subtitle = sprintf("%d observations; %d papers; %d site-years; top 10 papers = %s%% of observations",
                          nrow(L), nrow(per_paper), uniqueN(L$site_year), top_share(per_paper$n_obs, 10))) +
  theme_es + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
p_c <- heat(sg_region, "Region", "papers", "C. System group x region (unique papers)") + theme(legend.position = "none")
p4 <- (p_loc + labs(subtitle = NULL) | p_b) / (p_c | p_d) + plot_layout(heights = c(1.15, 1))
fig(p4, "08_evidence_sufficiency_4panel", 16, 13)

# ---- 12. Research-question feasibility ------------------------------------------------------
n_rep <- sum(cand$within_study_repeated); n_my <- sum(cand$multi_year); n_w <- sum(cand$water_contrast)
cl <- group_tab[system_group == "Cereal + legume"]
groups_5plus <- group_tab[papers >= 5 & !grepl("Unidentified", system_group), .N]
feas <- data.table(
  question = c("1. Does intercropping performance vary with climatic stress?",
               "2. Do functional intercropping systems respond differently to climatic stress?",
               "3. Are there climatic/environmental limits beyond which intercropping advantage declines?",
               "4. Does intercropping reduce or avoid productivity losses under stressful conditions?",
               "5. Can empirical intercropping relationships be transferred spatially across SSA?"),
  evidence_required = c(
    "Many independent site-years with year-specific LER spanning a wide range of stress, with enough site-years per stress level",
    "As (1), within EACH system group, and system groups observed under overlapping stress ranges and regions",
    "Dense evidence at the extremes of the stress gradient; enough independent site-years to estimate a threshold or slope change",
    "Paired observations of the same system under stressed and non-stressed conditions (multi-year, multi-season or water contrasts), ideally with sole-crop yields",
    "Broad, even geographic coverage; environmental conditions at prediction locations represented in the training data; out-of-region validation"),
  can_assess_now = c(
    sprintf("Number of independent site-years (%d) and how many are year-specific; geographic and temporal spread", uniqueN(L$site_year)),
    sprintf("Papers and site-years per system group; %d groups have 5+ papers; system-geography confounding", groups_5plus),
    "Only the volume and spread of independent evidence; not whether extremes are represented",
    sprintf("%d papers with a potential within-study repeated-condition comparison (%d multi-year, %d water/irrigation contrast)", n_rep, n_my, n_w),
    "Geographic coverage of papers and locations (maps, grid density, regional gaps)"),
  needs_climate_data = c(
    "Stress values per site-year; whether the stress gradient is wide and evenly sampled",
    "Whether system groups overlap in stress space, or system and climate are confounded",
    "Coverage of the stress extremes; whether extreme site-years come from few papers",
    "Whether repeated years/seasons actually differ in stress; which condition is the stressed one",
    "Environmental similarity of target areas to sampled sites (extrapolation risk)"),
  structural_assessment = c("Promising", "Uncertain", "Uncertain", "Uncertain", "Not currently supported"),
  explanation = c(
    sprintf("%d papers and %d site-years across %d countries give a reasonable base, but observations are clustered (top 10 papers = %s%% of observations) and some site-years pool several years, so the effective sample is the site-year count, not %d rows",
            nrow(per_paper), uniqueN(L$site_year), uniqueN(na.omit(L$Country)), top_share(per_paper$n_obs, 10), nrow(L)),
    sprintf("Evidence is dominated by cereal + legume (%d papers); most other groups have few papers and are geographically concentrated (see system x geography heatmap), so system and location effects may be hard to separate",
            cl$papers),
    "Thresholds need many independent site-years at the stress extremes; with evidence clustered in few papers and regions this is unlikely to be well identified, but it cannot be judged until stress is attached",
    sprintf("Feasible only for the %d papers with within-study repeated conditions; whether these contrasts involve real stress differences is unknown until climate is linked", n_rep),
    "Coverage is concentrated in a few countries and regions, with large gaps; spatial transfer would require extrapolation, which a structural count cannot justify"))
tab(feas, "09_research_question_feasibility")

# ---- 13. Assumptions and derived dataset -----------------------------------------------------
assumptions <- data.table(assumption = c(
  "Input dataset outputs/LER_all.csv is read only and not modified; every row is kept (missing fields exclude a row only from analyses that need that field).",
  "Missing country, coordinates, author and publication year for raw-only papers are filled from the ERA raw extraction (Data.Out) by paper x Site.ID; fills are listed in metadata_filled_from_raw.csv.",
  "Paper = ERA study code. Paper-site = Code x Site.ID as recorded (fertility levels or blocks recorded as sites count separately). Location = coordinates rounded to 0.01 degrees, pooled across papers.",
  "Site.ID values joined with '..' record several sites pooled into one observation; they are one paper-site and usually lack coordinates.",
  "Experiment year = calendar years in M.Year (four-digit years). A suffix such as .1/.2 or S/L is a season within the year. Two or more years in one M.Year (e.g. 2013.2014) is a multi-year average and counts as one site-year. M.Year without a calendar year (e.g. 'Year 1.Year 2') has no experiment year.",
  "Site-year = Code x Site.ID x M.Year as recorded; a multi-year average is one site-year.",
  "Crop combinations are order-independent, with spelling/case normalised (lookups/crop_type_lookup.csv). A single-crop product label is replaced by the Diversity combination when that names 2+ crops.",
  "System groups: any tree/perennial crop -> Tree-based; one crop -> Unidentified; 3+ crops -> Multi-crop; otherwise by the two crop types (lookups/system_group_crosswalk.csv).",
  "Regions are UN sub-regions (Natural Earth).",
  "Treatment information (practice codes, irrigation) is used only for rows where it describes the LER observation (compiled reported rows and rows calculated from yields). For reported_raw rows these fields were copied from another row of the paper, so only treatment text is used, and text cues are reported separately as unverified.",
  "Water practices = ERA practice codes under Blue or Green Water Management. Water/irrigation contrast = irrigated and rainfed observations, an irrigation difference within a comparison, a water practice in plist, or differing water practices between observations.",
  "Multi-season = two or more M.Year values within the same calendar year for one paper.",
  "No drought or climatic exposure is inferred; no environmental covariates are used."))
tab(assumptions, "assumptions")

derived <- L[, .(obs_id, paper, Author, Date, Country, ISO.3166.1.alpha.3, Region, Site.ID, paper_site, site_aggregated, Latitude, Longitude, Buffer,
                 location, M.Year, year_precision, years = sapply(years, paste, collapse = ";"), year_start, decade, site_year,
                 Product.Simple, Diversity, combination, combination_source, crop_types, system_group, T.Descrip, LER = MeanT, LER.Source,
                 Irrigation.C, Irrigation.T, water_practices_T = t_water, water_practice_in_plist = plist_water, water_text_cue,
                 treatment_fields_observation_specific = obs_specific_treatment)]
fwrite(derived, file.path(es_dir, "data", "LER_evidence_structure.csv"))
stopifnot(nrow(derived) == n_input)
