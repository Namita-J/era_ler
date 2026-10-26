# LER from ERA yield data: matching intercrop components to their monocrop references
#
# ERA stores pairwise comparisons (control arm vs treatment arm). Two facts make a purely
# row-based OR a purely treatment-based approach fail:
#  * an intercrop component is often only compared with a *base* monocrop (e.g. "MC PM F1"
#    vs "MM F0"), while the matching monocrop ("MM F1") only appears in another comparison;
#  * some papers reuse one control ID for several sole crops with different yields
#    (e.g. sole maize at three densities all coded "T0"), so only the comparison row says
#    which sole yield belongs to which intercrop.
# So for each intercrop component we (1) identify the monocrop treatment(s) of the same crop
# with the same experimental context, (2) use the yields paired in a direct ERA comparison
# row where one exists, and otherwise (3) look the monocrop yield up from the treatment table.
# Anything that cannot be resolved to exactly one (intercrop yield, monocrop yield) pair is
# flagged with a reason instead of being guessed.
#
# Requires: data.table, ic_codes (intercropping practice codes) in the calling environment.

library(data.table)

sole_kw <- "\\b(mono|monocrop|monocropping|monoculture|sole|pure|alone)\\b"
ic_kw   <- "inter[- ]?cr?op|\\bint\\b|\\binter\\b|\\bmixed\\b|\\bmixture\\b"   # incl. typos such as "intercop"
unit_kg <- c("kg/ha" = 1, "Mg/ha" = 1000, "t/ha" = 1000)
label_stop <- c("sole", "mono", "monocrop", "monocropping", "monoculture", "pure", "alone", "intercrop",
                "intercropped", "intercropping", "intercop", "int", "inter", "crop", "cropping", "with", "and", "the",
                "of", "plus", "row", "rows", "ha", "kg", "control")

norm_txt <- function(x) gsub("[^a-z0-9]", "", tolower(x))
sys_of   <- function(lbl) ifelse(grepl(">>", lbl, fixed = TRUE), sub("^.*?>>", "", lbl), lbl)
comp_of  <- function(lbl) ifelse(grepl(">>", lbl, fixed = TRUE), sub(">>.*$", "", lbl), lbl)
sym_diff <- function(a, b) union(setdiff(a, b), setdiff(b, a))

sys_parts <- function(s) {
  if (grepl("***", s, fixed = TRUE)) return(trimws(strsplit(s, "***", fixed = TRUE)[[1]]))
  if (grepl("..", s, fixed = TRUE) && !grepl("...", s, fixed = TRUE)) return(trimws(strsplit(s, "..", fixed = TRUE)[[1]]))
  trimws(s)
}

# Label tokens: letter+number codes (F0, N3, T1), numbers (10, 22.5) and words (low, sprayed).
# "CCF2" is split into "cc" + "f2" so it compares with "MC Cp F2".
mod_tokens <- function(s) {
  s <- gsub("([a-z])([a-z][0-9])", "\\1 \\2", tolower(s))
  tk <- regmatches(s, gregexpr("[a-z]?[0-9]+(\\.[0-9]+)?|[a-z]{2,}", s))[[1]]
  unique(setdiff(tk, label_stop))
}

code_set <- function(dt, cols) {
  m <- as.matrix(dt[, ..cols])
  apply(m, 1, function(r) paste(sort(unique(r[!is.na(r) & r != ""])), collapse = " "))
}

variety_ok <- function(a, b) {
  if (is.na(a) || is.na(b) || !nzchar(a) || !nzchar(b)) return(TRUE)
  a <- norm_txt(a); b <- norm_txt(b)
  grepl(a, b, fixed = TRUE) || grepl(b, a, fixed = TRUE)
}

# ---- 0. Outcome identity from the raw extraction -----------------------------------------
# The compiled table maps raw outcomes "Crop Yield..t/ha" and "Crop Yield..t/ha..LER" to the same
# Outcode (101). Returns the compiled row Index of yield rows whose raw outcome is LER-coded,
# matched on paper x data location x value.
ler_coded_index <- function(E, data_out) {
  normloc <- function(x) gsub("[[:space:]]", "", sub("^figure", "fig", sub("^table", "tab", tolower(x))))
  raw <- as.data.table(data_out)[grepl("[.][.]LER$", Out.Code.Joined) & !grepl("^Land Equivalent", Out.Code.Joined)]
  keys <- unique(raw[, paste(B.Code, normloc(ED.Data.Loc), round(as.numeric(ED.Mean.T), 3))])
  E[Outcode == "101" & paste(Code, normloc(DataLoc), round(MeanT, 3)) %in% keys, Index]
}

# ---- 1. Comparison rows and treatment table ----------------------------------------------
build_tables <- function(E, codes, exclude_index = integer(0)) {
  y <- E[Code %in% codes & Outcode == "101" & !grepl("-", Product.Simple)]
  # rows whose raw ERA outcome is in LER units (e.g. "Crop Yield..t/ha..LER": reported LERs or
  # pLERs entered under the yield outcome) are not yields; they are logged and not used
  ler_coded_log <- y[Index %in% exclude_index, .(Index, Code, Site.ID, M.Year, Crop = Product.Simple, TID, T.Descrip, MeanT, DataLoc)]
  y <- y[!Index %in% exclude_index]
  y[, Ccodes := code_set(y, paste0("C", 1:13))]
  y[, Tcodes := code_set(y, paste0("T", 1:13))]
  rows <- y[, .(Index, Code, Site.ID, M.Year, Crop = Product.Simple, Units, Variety,
                CID, TID, MeanC, MeanT, DataLoc)]

  arms <- rbind(
    y[, .(Index, Code, Site.ID, M.Year, Crop = Product.Simple, Units, Variety,
          id = CID, label = C.Descrip, yield = MeanC, codes = Ccodes)],
    y[, .(Index, Code, Site.ID, M.Year, Crop = Product.Simple, Units, Variety,
          id = TID, label = T.Descrip, yield = MeanT, codes = Tcodes)]
  )
  arms <- arms[!is.na(yield) & !is.na(id)]
  # a 0 recorded for a treatment that also has non-zero yields in the same block is a
  # placeholder (e.g. an empty control cell), not an observation
  arms[, zero_placeholder := yield == 0 & any(yield > 0), by = .(Code, Site.ID, M.Year, Crop, id)]
  zero_log <- unique(arms[zero_placeholder == TRUE, .(Code, Site.ID, M.Year, Crop, id, label)])
  arms <- arms[zero_placeholder == FALSE]

  tr <- arms[, .(
    label    = label[1],
    codes    = codes[1],
    n_codes  = uniqueN(codes),
    Variety  = paste(unique(na.omit(Variety)), collapse = "; "),
    Units    = paste(unique(Units), collapse = "; "),
    yields   = list(sort(unique(round(yield, 6)))),
    n_yields = uniqueN(round(yield, 6)),
    rows     = paste(unique(Index), collapse = ",")
  ), by = .(Code, Site.ID, M.Year, Crop, id)]
  tr[, yield := sapply(yields, function(v) if (length(v) == 1) v else NA_real_)]
  tr[, `:=`(comp = comp_of(label), sys = sys_of(label))]
  tr[, code_vec := lapply(strsplit(codes, " "), function(v) v[nzchar(v)])]
  tr[, has_ic_code := sapply(code_vec, function(v) any(v %in% ic_codes))]
  tr[, `:=`(lab_sole = grepl(sole_kw, comp, ignore.case = TRUE),
            lab_ic   = grepl(ic_kw,   comp, ignore.case = TRUE))]
  # a label shared by treatments of 2+ crops is an intercrop system, unless it is labelled sole
  # (sole crops of different species can share a generic label such as "sole unsprayed")
  tr[, sys_multi := uniqueN(Crop) > 1, by = .(Code, Site.ID, M.Year, sys)]
  tr[, is_ic := has_ic_code | (sys_multi & !lab_sole)]
  # practice codes and label disagree on whether this arm is a monocrop or an intercrop
  tr[, label_code_conflict := (is_ic & lab_sole & !lab_ic) | (!is_ic & lab_ic & !lab_sole)]
  # a monocrop-coded treatment whose label is identical to an intercrop component label of the
  # same crop (e.g. control "MZ- 2:2 rows-maize" = maize component of the 2:2 intercrops) is an
  # intercrop component used as a reference, not a sole crop
  tr[, ic_comp_labels := list(list(unique(norm_txt(comp[is_ic & !lab_sole])))), by = .(Code, Site.ID, M.Year, Crop)]
  tr[, label_is_ic_component := !is_ic & mapply(function(c, l) norm_txt(c) %in% l, comp, ic_comp_labels)]
  tr[label_is_ic_component == TRUE, label_code_conflict := TRUE]
  tr[, ic_comp_labels := NULL]
  tr[, ctx_codes := lapply(code_vec, function(v) setdiff(v, ic_codes))]
  tr[, full_tokens := lapply(label, mod_tokens)]
  list(rows = rows, treatments = tr, zero_log = zero_log, ler_coded_log = ler_coded_log)
}

# ---- 2. Match one intercrop component ---------------------------------------------------
# ic: one intercrop treatment; monos: valid monocrop treatments of the same crop in the same
# paper x site x year; direct: comparison rows linking ic to one of those monos, oriented as
# (ic_y, mono_id, mono_y); factors: practice codes that vary between treatments of this crop.
match_component <- function(ic, monos, direct, factors) {
  out <- list(mono_id = NA_character_, mono_label = NA_character_, mono_codes = NA_character_,
              ic_y = NA_real_, mono_y = NA_real_, match_basis = NA_character_,
              factor_diff = NA_character_, design_diff = NA_character_,
              reason = NA_character_, candidates = NA_character_)
  if (nrow(monos) == 0) { out$reason <- "no_matching_monocrop"; return(out) }
  ic_ctx <- ic$ctx_codes[[1]]; ic_tok <- ic$full_tokens[[1]]; ic_dist <- ic$ic_dist[[1]]
  monos <- copy(monos)[sapply(Variety, variety_ok, b = ic$Variety)]
  if (nrow(monos) == 0) { out$reason <- "treatment_mismatch"; return(out) }

  # label modifiers that distinguish the candidate monocrops from each other (F0/F1, low/high,
  # sprayed/unsprayed, Soybean-10/-20); tokens shared by all candidates carry no information.
  # Two-way: the monocrop's modifiers must appear in the intercrop label, and the intercrop's own
  # modifiers (e.g. "mulch" in "CT + Mulch Cowpea intercrop") must appear in the monocrop label.
  common <- Reduce(intersect, monos$full_tokens)
  monos[, label_pos := sapply(full_tokens, function(t) {
    d <- setdiff(t, common)
    length(d) > 0 && all(d %in% ic_tok) && all(ic_dist %in% t)
  })]
  monos[, f_diff := lapply(ctx_codes, function(v) intersect(sym_diff(v, ic_ctx), factors))]
  monos[, n_f_diff := lengths(f_diff)]

  # same experimental context: label evidence first (it can override inconsistent practice
  # coding), otherwise agreement on the experimental-factor practice codes
  cand <- if (any(monos$label_pos)) monos[label_pos == TRUE] else monos[n_f_diff == 0]
  if (nrow(cand) == 0) { out$reason <- "treatment_mismatch"; return(out) }
  cand <- cand[n_f_diff == min(n_f_diff)]

  d <- unique(direct[mono_id %in% cand$id, .(ic_y = round(ic_y, 6), mono_y = round(mono_y, 6), mono_id)])
  d <- d[mono_y > 0]
  if (nrow(d) > 0) {
    out$candidates <- paste0(d$mono_id, ":", d$ic_y, "/", d$mono_y, collapse = "; ")
    # an intercrop treatment with several yields in one paper x site x year (unrecorded
    # seasons/sites/harvests) cannot be paired safely with the other component crop
    if (ic$n_yields > 1) { out$reason <- if (nrow(d) > 1) "multiple_yield_pairs" else "multiple_intercrop_yields"; return(out) }
    if (nrow(d) > 1) {
      out$reason <- if (uniqueN(d$mono_id) > 1) "multiple_matching_monocrops" else "multiple_yield_pairs"
      return(out)
    }
    sel <- cand[id == d$mono_id]
    out$ic_y <- d$ic_y; out$mono_y <- d$mono_y; out$match_basis <- "direct_comparison"
  } else {
    out$candidates <- paste0(cand$id, "=", sapply(cand$yields, paste, collapse = "/"), collapse = "; ")
    if (nrow(cand) > 1)      { out$reason <- "multiple_matching_monocrops"; return(out) }
    if (cand$n_yields > 1)   { out$reason <- "multiple_matching_monocrops"; return(out) }
    if (ic$n_yields > 1)     { out$reason <- "multiple_intercrop_yields";   return(out) }
    sel <- cand
    out$ic_y <- ic$yield; out$mono_y <- cand$yield; out$match_basis <- "treatment_lookup"
  }
  out$mono_id <- sel$id; out$mono_label <- sel$label; out$mono_codes <- sel$codes
  out$factor_diff <- paste(sel$f_diff[[1]], collapse = " ")
  out$design_diff <- paste(setdiff(sym_diff(sel$ctx_codes[[1]], ic_ctx), factors), collapse = " ")
  if (sel$label_pos) out$match_basis <- paste0(out$match_basis, "+label")
  out
}

# ---- 3. Assemble intercrop systems and compute LER -----------------------------------
calc_ler_systems <- function(E, codes, exclude_index = integer(0)) {
  tb <- build_tables(E, codes, exclude_index)
  tr <- tb$treatments; rows <- tb$rows
  crop_words <- unique(tolower(unlist(strsplit(unique(E$Product.Simple), "[^A-Za-z]+"))))
  crop_words <- crop_words[nchar(crop_words) >= 4]

  ic_tr <- tr[is_ic == TRUE]
  ok_mono <- tr[is_ic == FALSE & label_code_conflict == FALSE]
  # an intercrop component's own modifiers: tokens of its component label (before ">>", so the
  # partner crop's variety is not included) that differ between intercrop components of that crop
  ic_tr[, comp_tok := lapply(comp, mod_tokens)]
  ic_tr[, ic_dist := { cm <- Reduce(intersect, comp_tok); list(lapply(comp_tok, function(t) setdiff(t, cm))) },
        by = .(Code, Site.ID, M.Year, Crop)]
  # ...restricted to the sole-crop vocabulary of that crop (words used in any monocrop-coded label,
  # including mislabelled ones): only modifiers by which sole crops are distinguished must match
  vocab <- tr[is_ic == FALSE, .(vocab = list(unique(unlist(full_tokens)))), by = .(Code, Site.ID, M.Year, Crop)]
  ic_tr <- merge(ic_tr, vocab, by = c("Code", "Site.ID", "M.Year", "Crop"), all.x = TRUE)
  ic_tr[, ic_dist := mapply(function(d, v) intersect(d, v), ic_dist, vocab, SIMPLIFY = FALSE)]

  # comparison rows oriented as intercrop arm vs valid monocrop arm (either direction)
  key <- c("Code", "Site.ID", "M.Year", "Crop")
  r1 <- merge(rows, ic_tr[, .(Code, Site.ID, M.Year, Crop, TID = id)], by = c(key, "TID"))
  r1 <- merge(r1, ok_mono[, .(Code, Site.ID, M.Year, Crop, CID = id)], by = c(key, "CID"))
  r2 <- merge(rows, ic_tr[, .(Code, Site.ID, M.Year, Crop, CID = id)], by = c(key, "CID"))
  r2 <- merge(r2, ok_mono[, .(Code, Site.ID, M.Year, Crop, TID = id)], by = c(key, "TID"))
  direct <- rbind(r1[, .(Code, Site.ID, M.Year, Crop, ic_id = TID, mono_id = CID, ic_y = MeanT, mono_y = MeanC)],
                  r2[, .(Code, Site.ID, M.Year, Crop, ic_id = CID, mono_id = TID, ic_y = MeanC, mono_y = MeanT)])
  direct <- direct[!is.na(ic_y) & !is.na(mono_y)]

  # experimental factors per paper x site x year x crop: practice codes that vary between
  # intercrop treatments, or between monocrop treatments, of that crop
  varying <- function(sets) if (length(sets) < 2) character(0) else setdiff(Reduce(union, sets), Reduce(intersect, sets))
  fac <- tr[label_code_conflict == FALSE,
            .(factors = list(union(varying(ctx_codes[is_ic]), varying(ctx_codes[!is_ic])))), by = key]

  # components expected in each system ID: parts that are an intercrop component label
  # somewhere in the paper, or that name a crop (excludes tags such as "Zizi" for a tree)
  paper_comps <- ic_tr[, .(comps = list(unique(norm_txt(comp)))), by = Code]
  sys_tab <- merge(unique(ic_tr[, .(Code, Site.ID, M.Year, sys)]), paper_comps, by = "Code")
  sys_tab[, n_expected := mapply(function(s, cps) {
    p <- sys_parts(s)
    sum(norm_txt(p) %in% cps | sapply(tolower(p), function(x) any(sapply(crop_words, grepl, x = x, fixed = TRUE))))
  }, sys, comps)]
  sys_tab[, comps := NULL]

  comp <- rbindlist(lapply(seq_len(nrow(ic_tr)), function(i) {
    ic <- ic_tr[i]
    blk <- list(ic$Code, ic$Site.ID, ic$M.Year, ic$Crop)
    monos <- ok_mono[Code == blk[[1]] & Site.ID == blk[[2]] & M.Year == blk[[3]] & Crop == blk[[4]]]
    dr <- direct[Code == blk[[1]] & Site.ID == blk[[2]] & M.Year == blk[[3]] & Crop == blk[[4]] & ic_id == ic$id]
    f <- fac[Code == blk[[1]] & Site.ID == blk[[2]] & M.Year == blk[[3]] & Crop == blk[[4]], factors]
    m <- if (ic$label_code_conflict) list(reason = "label_code_conflict")
         else match_component(ic, monos, dr, if (length(f)) f[[1]] else character(0))
    as.data.table(c(list(Code = ic$Code, Site.ID = ic$Site.ID, M.Year = ic$M.Year, sys = ic$sys,
                         Crop = ic$Crop, ic_id = ic$id, ic_label = ic$label, ic_codes = ic$codes,
                         ic_units = ic$Units, ic_rows = ic$rows,
                         ic_yields_all = paste(ic$yields[[1]], collapse = "/")), m))
  }), fill = TRUE)

  # pLER at full precision (monocrop and intercrop yields come from the same row or the same
  # unit-consistent treatment table, so units are shared)
  comp[, pLER := ic_y / mono_y]
  comp[is.na(reason) & (is.na(mono_y) | mono_y <= 0), reason := "zero_monocrop_yield"]
  comp[is.na(reason) & ic_y < 0, reason := "negative_intercrop_yield"]   # e.g. a digitising artefact below zero
  comp[, n_ic_crop := uniqueN(ic_id), by = .(Code, Site.ID, M.Year, sys, Crop)]
  comp[is.na(reason) & n_ic_crop > 1, reason := "multiple_intercrop_treatments"]

  sys <- comp[, .(
    n_crops   = uniqueN(Crop),
    crops     = paste(sort(unique(Crop)), collapse = "-"),
    n_bad     = sum(!is.na(reason)),
    reasons   = paste(sort(unique(na.omit(reason))), collapse = "; "),
    pLER      = paste0(Crop, "=", sprintf("%.3f", pLER), collapse = "; "),
    LER       = sum(pLER),
    label_override = any(grepl("label", match_basis) & nzchar(factor_diff), na.rm = TRUE),
    design_diff    = any(nzchar(design_diff), na.rm = TRUE),
    ic_ids    = paste(ic_id, collapse = "-"),
    mono_ids  = paste(mono_id, collapse = "-")
  ), by = .(Code, Site.ID, M.Year, sys)]
  sys <- merge(sys, sys_tab, by = c("Code", "Site.ID", "M.Year", "sys"), all.x = TRUE)
  sys[, status := fcase(
    n_bad > 0,                                    reasons,
    n_crops < 2,                                  "missing_intercrop_component",
    !is.na(n_expected) & n_crops < n_expected,    "missing_intercrop_component",
    default = "OK")]
  sys[status != "OK", LER := NA_real_]
  # QC flags on accepted systems (not reasons for rejection):
  #  factor_codes_differ_label_match: treatment labels match but ERA practice codes for an experimental factor disagree
  #  design_management_differs: sole crop and intercrop differ only in management that is constant across the experiment
  sys[, qc_flag := fifelse(status != "OK", NA_character_,
                   fifelse(label_override & design_diff, "factor_codes_differ_label_match; design_management_differs",
                   fifelse(label_override, "factor_codes_differ_label_match",
                   fifelse(design_diff, "design_management_differs", NA_character_))))]
  list(components = comp, systems = sys, treatments = tr, direct = direct, zero_log = tb$zero_log,
       ler_coded_log = tb$ler_coded_log)
}

# ---- Reported LER from the raw extraction ---------------------------------------------
# The compiled table only keeps LER (Outcode 103) where it was entered as a comparison. Raw
# Data.Out also holds LER recorded as an outcome in its own right ("Land Equivalent Ratio"),
# and LER entered under the yield outcome with "LER" in Species/Item/Group (Out.Group).
# Returns one row per reported LER observation, with partial (per-crop) LERs separated out.
reported_ler_from_raw <- function(data_out, out_out) {
  O <- as.data.table(out_out)
  ler_out <- O[(Out.Subind == "Land Equivalent Ratio" & !grepl("partial", fifelse(is.na(Out.Group), "", Out.Group), ignore.case = TRUE)) |
               grepl("^LER$", fifelse(is.na(Out.Group), "", Out.Group)),
               .(B.Code, Out.Code.Joined, LER.Group = Out.Group)]
  partial_out <- O[Out.Subind == "Land Equivalent Ratio" & grepl("partial", fifelse(is.na(Out.Group), "", Out.Group), ignore.case = TRUE),
                   .(B.Code, Out.Code.Joined)]
  R <- merge(as.data.table(data_out), unique(ler_out), by = c("B.Code", "Out.Code.Joined"))
  R[, LER := suppressWarnings(as.numeric(ED.Mean.T))]
  R <- R[!is.na(LER)]
  R[, IN := fifelse(is.na(IN.Level.Name), "", IN.Level.Name)]
  # identity of an observation: paper x site x time x treatment x system x all treatment level names
  lvl <- setdiff(grep("[.]Level[.]Name$", names(R), value = TRUE), c("IN.Level.Name", "R.Level.Name"))
  R[, levels_key := do.call(paste, c(lapply(.SD, function(v) fifelse(is.na(v), "", as.character(v))), sep = "|")), .SDcols = lvl]
  # partial LERs: 2+ single-crop products of one system carrying different values
  R[, single_crop := !is.na(P.Product) & !grepl("-", P.Product)]
  R[, partial := sum(single_crop) >= 2 & uniqueN(P.Product[single_crop]) >= 2 & uniqueN(round(LER[single_crop], 4)) >= 2,
    by = .(B.Code, Site.ID, Time, IN, levels_key, Out.Code.Joined)]
  partial_log <- rbind(R[partial == TRUE, .(B.Code, Site.ID, Time, T.Name, IN.Level.Name, P.Product, LER, Out.Code.Joined, ED.Data.Loc, why = "per-crop values (partial LER)")],
                       merge(as.data.table(data_out), partial_out, by = c("B.Code", "Out.Code.Joined"))[
                         , .(B.Code, Site.ID, Time, T.Name, IN.Level.Name, P.Product, LER = suppressWarnings(as.numeric(ED.Mean.T)), Out.Code.Joined, ED.Data.Loc, why = "outcome group 'Partial'")])
  R <- R[partial == FALSE]
  # the same system LER entered once per component crop is one observation
  R <- R[!duplicated(R[, .(B.Code, Site.ID, Time, IN, T.Name, levels_key, Out.Code.Joined, round(LER, 6))])]
  R[, n_values_system := uniqueN(round(LER, 4)), by = .(B.Code, Site.ID, Time, IN, T.Name, levels_key, Out.Code.Joined)]
  out <- R[, .(Code = B.Code, Site.ID, M.Year = gsub("..", ".", Time, fixed = TRUE), Time, T.Name, IN.Level.Name,
               Product = fifelse(is.na(P.Product), IN.Prods.All, P.Product), LER, Out.Code.Joined, LER.Group, DataLoc = ED.Data.Loc,
               levels_key, n_values_system)]
  list(reported = out, partial_log = partial_log)
}
