# PROTOTYPE: LER from the raw ERA extraction (Data.Out) instead of the compiled comparisons.
#
# Raw Data.Out holds one record per treatment observation, with the component treatment
# (T.Name), its intercrop system (IN.Level.Name), practice codes (T.Codes, IN.T.Codes), the
# original outcome identity (Out.Code.Joined) and all treatment level names. This keeps
# records and identities that the compiled comparison table drops or distorts.
#
# Reuses the matching rules in ler_matching.R (match_component, mod_tokens, ...).
# Requires: data.table, ic_codes, source("r_scripts/ler_matching.R").

library(data.table)

split_codes <- function(x) { v <- unlist(strsplit(ifelse(is.na(x), "", x), "-", fixed = TRUE)); v[nzchar(v)] }

# ---- 1. Treatment table from raw records ------------------------------------------------
build_raw_treatments <- function(data_out, codes) {
  D <- as.data.table(data_out)[B.Code %in% codes]
  # yield outcome only; LER-coded outcomes ("..LER") are kept separately as reported LER
  ler_reported <- D[grepl("[.][.]LER$", Out.Code.Joined)]
  D <- D[Out.Subind == "Crop Yield" & !grepl("[.][.]LER$", Out.Code.Joined) & !is.na(P.Product)]
  D[, yield := suppressWarnings(as.numeric(ED.Mean.T))]
  D <- D[!is.na(yield)]
  D[, IN.Level.Name := fifelse(is.na(IN.Level.Name), "", IN.Level.Name)]
  # a 0 recorded for a treatment that also has non-zero yields is a placeholder (as in ler_matching.R)
  D[, zero_placeholder := yield == 0 & any(yield > 0), by = .(B.Code, Site.ID, Time, P.Product, T.Name, IN.Level.Name)]
  D <- D[zero_placeholder == FALSE]

  # treatment identity = treatment name + system + every treatment level name / variety field,
  # so treatments that share a name but differ in another recorded level stay distinct
  lvl_cols <- intersect(c(grep("\\.Level\\.Name$", names(D), value = TRUE), "V.Var", "Structure.Comb"), names(D))
  lvl_cols <- setdiff(lvl_cols, c("IN.Level.Name", "R.Level.Name"))
  D[, levels_key := do.call(paste, c(lapply(.SD, function(v) fifelse(is.na(v), "", as.character(v))), sep = "|")), .SDcols = lvl_cols]

  tr <- D[, .(
    codes    = paste(sort(unique(c(split_codes(T.Codes[1]), split_codes(IN.T.Codes[1])))), collapse = " "),
    Variety  = paste(unique(na.omit(c(V.Var, V.Level.Name))), collapse = "; "),
    Units    = paste(unique(Out.Unit), collapse = "; "),
    yields   = list(sort(unique(round(yield, 6)))),
    n_yields = uniqueN(round(yield, 6)),
    data_loc = paste(unique(ED.Data.Loc), collapse = "; "),
    n_records = .N
  ), by = .(Code = B.Code, Site.ID, M.Year = Time, Crop = P.Product, T.Name, IN.Level.Name, levels_key)]
  tr[, yield := sapply(yields, function(v) if (length(v) == 1) v else NA_real_)]
  tr[, id := paste0("R", .I)]
  tr[, `:=`(comp = T.Name, sys = fifelse(nzchar(IN.Level.Name), IN.Level.Name, T.Name),
            label = fifelse(nzchar(IN.Level.Name), paste0(T.Name, ">>", IN.Level.Name), T.Name))]
  tr[, code_vec := lapply(strsplit(codes, " "), function(v) v[nzchar(v)])]
  tr[, ctx_codes := lapply(code_vec, function(v) setdiff(v, ic_codes))]

  # an intercrop component: its system holds 2+ crops at this site and time
  tr[, n_sys_crops := uniqueN(Crop), by = .(Code, Site.ID, M.Year, IN.Level.Name)]
  tr[, is_ic := nzchar(IN.Level.Name) & n_sys_crops >= 2]
  tr[, `:=`(lab_sole = grepl(sole_kw, comp, ignore.case = TRUE), lab_ic = grepl(ic_kw, comp, ignore.case = TRUE))]
  # the component's treatment name must be one of the parts of its system (catches records
  # filed under the wrong system, e.g. a "within row" maize filed under a 1:1 system)
  tr[, component_system_mismatch := is_ic & !mapply(function(t, s) norm_txt(t) %in% norm_txt(sys_parts(s)), T.Name, IN.Level.Name)]
  tr[, label_code_conflict := (is_ic & lab_sole & !lab_ic) | component_system_mismatch]
  # safeguard (as in ler_matching.R): a record outside any intercrop system whose treatment name is
  # an intercrop component name of the same crop (e.g. "MZ- 2:2 rows-maize") is not a sole crop
  tr[, ic_comp_labels := list(list(unique(norm_txt(comp[is_ic & !lab_sole])))), by = .(Code, Site.ID, M.Year, Crop)]
  tr[, label_is_ic_component := !is_ic & mapply(function(c, l) norm_txt(c) %in% l, comp, ic_comp_labels)]
  tr[label_is_ic_component == TRUE, label_code_conflict := TRUE]
  tr[, ic_comp_labels := NULL]
  tr[, full_tokens := lapply(label, mod_tokens)]
  list(treatments = tr, ler_reported = ler_reported)
}

# ---- 2. Systems and LER -------------------------------------------------------------
calc_ler_raw <- function(data_out, codes) {
  tb <- build_raw_treatments(data_out, codes)
  tr <- tb$treatments
  ic_tr <- tr[is_ic == TRUE]
  ok_mono <- tr[is_ic == FALSE & label_code_conflict == FALSE]
  if (nrow(ic_tr) == 0) return(list(components = data.table(), systems = data.table(), treatments = tr, ler_reported = tb$ler_reported))

  ic_tr[, comp_tok := lapply(comp, mod_tokens)]
  ic_tr[, ic_dist := { cm <- Reduce(intersect, comp_tok); list(lapply(comp_tok, function(t) setdiff(t, cm))) },
        by = .(Code, Site.ID, M.Year, Crop)]
  vocab <- tr[is_ic == FALSE, .(vocab = list(unique(unlist(full_tokens)))), by = .(Code, Site.ID, M.Year, Crop)]
  ic_tr <- merge(ic_tr, vocab, by = c("Code", "Site.ID", "M.Year", "Crop"), all.x = TRUE)
  ic_tr[, ic_dist := mapply(function(d, v) intersect(d, v), ic_dist, vocab, SIMPLIFY = FALSE)]

  key <- c("Code", "Site.ID", "M.Year", "Crop")
  varying <- function(sets) if (length(sets) < 2) character(0) else setdiff(Reduce(union, sets), Reduce(intersect, sets))
  fac <- tr[label_code_conflict == FALSE,
            .(factors = list(union(varying(ctx_codes[is_ic]), varying(ctx_codes[!is_ic])))), by = key]
  no_direct <- data.table(mono_id = character(0), ic_y = numeric(0), mono_y = numeric(0))

  comp <- rbindlist(lapply(seq_len(nrow(ic_tr)), function(i) {
    ic <- ic_tr[i]
    monos <- ok_mono[Code == ic$Code & Site.ID == ic$Site.ID & M.Year == ic$M.Year & Crop == ic$Crop]
    f <- fac[Code == ic$Code & Site.ID == ic$Site.ID & M.Year == ic$M.Year & Crop == ic$Crop, factors]
    m <- if (ic$component_system_mismatch) list(reason = "component_system_mismatch")
         else if (ic$label_code_conflict) list(reason = "label_code_conflict")
         else match_component(ic, monos, no_direct, if (length(f)) f[[1]] else character(0))
    as.data.table(c(list(Code = ic$Code, Site.ID = ic$Site.ID, M.Year = ic$M.Year, sys = ic$IN.Level.Name,
                         Crop = ic$Crop, ic_id = ic$id, ic_label = ic$label, ic_codes = ic$codes,
                         ic_units = ic$Units, ic_yields_all = paste(ic$yields[[1]], collapse = "/")), m))
  }), fill = TRUE)

  comp[, pLER := ic_y / mono_y]
  comp[is.na(reason) & (is.na(mono_y) | mono_y <= 0), reason := "zero_monocrop_yield"]
  comp[is.na(reason) & ic_y < 0, reason := "negative_intercrop_yield"]   # e.g. a digitising artefact below zero
  comp[, n_ic_crop := uniqueN(ic_id), by = .(Code, Site.ID, M.Year, sys, Crop)]
  comp[is.na(reason) & n_ic_crop > 1, reason := "multiple_intercrop_treatments"]

  crop_words <- unique(tolower(unlist(strsplit(unique(tr$Crop), "[^A-Za-z]+")))); crop_words <- crop_words[nchar(crop_words) >= 4]
  paper_comps <- ic_tr[, .(comps = list(unique(norm_txt(comp)))), by = Code]
  sys_tab <- merge(unique(ic_tr[, .(Code, Site.ID, M.Year, sys = IN.Level.Name)]), paper_comps, by = "Code")
  sys_tab[, n_expected := mapply(function(s, cps) {
    p <- sys_parts(s); sum(norm_txt(p) %in% cps | sapply(tolower(p), function(x) any(sapply(crop_words, grepl, x = x, fixed = TRUE))))
  }, sys, comps)][, comps := NULL]

  sys <- comp[, .(
    n_crops = uniqueN(Crop), crops = paste(sort(unique(Crop)), collapse = "-"),
    n_bad = sum(!is.na(reason)), reasons = paste(sort(unique(na.omit(reason))), collapse = "; "),
    pLER = paste0(Crop, "=", sprintf("%.3f", pLER), collapse = "; "), LER = sum(pLER),
    label_override = any(grepl("label", match_basis) & nzchar(factor_diff), na.rm = TRUE),
    design_diff = any(nzchar(design_diff), na.rm = TRUE),
    ic_ids = paste(ic_id, collapse = "-"), mono_ids = paste(mono_id, collapse = "-")
  ), by = .(Code, Site.ID, M.Year, sys)]
  sys <- merge(sys, sys_tab, by = c("Code", "Site.ID", "M.Year", "sys"), all.x = TRUE)
  sys[, status := fcase(n_bad > 0, reasons, n_crops < 2, "missing_intercrop_component",
                        !is.na(n_expected) & n_crops < n_expected, "missing_intercrop_component", default = "OK")]
  sys[status != "OK", LER := NA_real_]
  sys[, qc_flag := fifelse(status != "OK", NA_character_,
                   fifelse(label_override & design_diff, "factor_codes_differ_label_match; design_management_differs",
                   fifelse(label_override, "factor_codes_differ_label_match",
                   fifelse(design_diff, "design_management_differs", NA_character_))))]
  list(components = comp, systems = sys, treatments = tr, ler_reported = tb$ler_reported)
}
