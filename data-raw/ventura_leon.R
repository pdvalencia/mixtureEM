# Prepare the bundled `ventura_leon` dataset from Ventura-Leon et al. (2025),
# "Exploring Infidelity Behavior Patterns in a Sample of Peruvian Young
# Adults: A Latent Class Analysis", Journal of Marital and Family Therapy,
# 51:e70066. https://doi.org/10.1111/jmft.70066
#
# Source: OSF repository https://osf.io/8csr9/ (data collected by the study
# authors; the OSF project itself has no license flag set, but Pablo D.
# Valencia -- a co-author of the paper and of this package -- has authorized
# bundling this dataset in mixtureEM on behalf of the author team).
#
# The raw file has 16 binary indicators (item1-item16) from the 20-item
# Multidimensional Infidelity Inventory-S (MII-S; Romero-Palencia et al.
# 2007) plus six sociodemographic variables. Recoded here to English
# variable names; values are otherwise used as-is (already binary/complete,
# no missing data).

raw <- readxl::read_excel("data-raw/ventura_leon/ventura_leon_raw.xlsx")

items <- as.data.frame(lapply(raw[paste0("item", 1:16)], as.integer))

# Informative indicator names, in the raw item1-item16 order.
item_names <- c(
  "flirting", "romantic_partners", "emotional_bond", "romantic_involvement",
  "loved_another", "in_love", "thoughts", "interest",
  "sexual_relations", "sexual_contact", "desired_relations", "desired_contact",
  "sexual_fantasies", "attraction", "had_sex", "desired_sex"
)
names(items) <- item_names

ventura_leon <- data.frame(
  sex = factor(raw$Sexo, levels = c("Varón", "Mujer"), labels = c("Male", "Female")),
  age = as.numeric(raw$Edad),
  birthplace = factor(raw$Lugar_nacimiento,
                       levels = c("Lima", "Fuera de Lima"),
                       labels = c("Lima", "Outside Lima")),
  sexual_orientation = factor(raw$Orientacion_sexual,
                               levels = c("Heterosexual", "No_heterosexual"),
                               labels = c("Heterosexual", "Not heterosexual")),
  relationship_duration = factor(raw$Tiempo_en_relacion,
                                  levels = c("Poco tiempo", "Mucho tiempo"),
                                  labels = c("Short", "Long")),
  relationship_type = factor(raw$Tipo_relacion,
                              levels = c("Casados", "Convivienetes", "Enamorados",
                                         "Novios", "Salientes"),
                              labels = c("Married", "Cohabiting", "In love",
                                         "Engaged", "Dating")),
  items
)

# Sanity checks against Table 1 (sociodemographics) and Table 3 (item total
# endorsement proportions) of Ventura-Leon et al. (2025).
stopifnot(nrow(ventura_leon) == 400)
stopifnot(all(table(ventura_leon$sex) == c(Male = 99, Female = 301)))
stopifnot(all(table(ventura_leon$birthplace) == c(Lima = 312, `Outside Lima` = 88)))
stopifnot(all(table(ventura_leon$sexual_orientation) ==
                c(Heterosexual = 348, `Not heterosexual` = 52)))
stopifnot(all(table(ventura_leon$relationship_duration) ==
                c(Short = 307, Long = 93)))
stopifnot(all(table(ventura_leon$relationship_type) ==
                c(Married = 36, Cohabiting = 34, `In love` = 222,
                  Engaged = 52, Dating = 56)))
stopifnot(abs(mean(ventura_leon$age) - 25.28) < 0.01)
stopifnot(abs(sd(ventura_leon$age) - 7.24) < 0.01)

item_totals <- round(sapply(ventura_leon[item_names], mean), 2)
published_totals <- c(0.50, 0.38, 0.44, 0.36, 0.31, 0.32, 0.58, 0.48, 0.20,
                       0.20, 0.28, 0.27, 0.24, 0.55, 0.16, 0.28)
stopifnot(all(item_totals == published_totals))

usethis::use_data(ventura_leon, overwrite = TRUE)
