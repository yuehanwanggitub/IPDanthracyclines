library(survival)
library(dplyr)
library(timereg)


# Base model
#-----------
model.base <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation), data = IPD.nN.n0, weights = weight)


# Selection procedure to evaluate other chemotherapeutic agents
#--------------------------------------------------------------
# Epipodophyllotoxins
model.Epipo <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + as.factor(Epipodophyllotoxins), data = IPD.nN.n0, weights = weight)

# Vinca alkaloids
model.Vinca <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + as.factor(VincaAlkaloids), data = IPD.nN.n0, weights = weight)

# Platinum compounds
model.Plat <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + as.factor(PlatinumCompounds), data = IPD.nN.n0, weights = weight)


# Antimetabolites
model.Antime <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + as.factor(Antimetabolites), data = IPD.nN.n0, weights = weight)


# Model 1
#--------
model.1 <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat) + as.factor(Daunorubicin.cat)
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat), data = IPD.nN.n0, weights = weight)
# Proportional hazards assumption
test.ph <- cox.zph(model.1)
test.ph


# Model 2
#--------
model.2 <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat), data = IPD.nN.n0.2D, weights = weight)
# Proportional hazards assumption
test.ph <- cox.zph(model.2)
test.ph


# Multiplicative interaction
#---------------------------
# Interaction: Doxorubicin dose * Chest radiotherapy (yes/no)
model.doxo.chestRT <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRT) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Doxorubicin.Doseall.per100:as.factor(ChestRT), data = IPD.nN.n0.2D, weights = weight)

# Interaction: Daunorubicin dose * Chest radiotherapy (yes/no)
model.dauno.chestRT <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRT) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Daunorubicin.Doseall.per100:as.factor(ChestRT), data = IPD.nN.n0.2D, weights = weight)

# Interaction: Doxorubicin dose * Chest radiotherapy field/dose combination
model.doxo.chestRTfielddose <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Doxorubicin.Doseall.per100:as.factor(ChestRTfield.totaldose.com), data = IPD.nN.n0.2D, weights = weight)

# Interaction: Daunorubicin dose * Chest radiotherapy field/dose combination
model.dauno.chestRTfielddose <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Daunorubicin.Doseall.per100:as.factor(ChestRTfield.totaldose.com), data = IPD.nN.n0.2D, weights = weight)

# Interaction: Doxorubicin dose * Age at childhood cancer diagnosis
model.doxo.ccage <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Doxorubicin.Doseall.per100:as.factor(ccage.cat), data = IPD.nN.n0.2D, weights = weight)

# Interaction: Daunorubicin dose * Age at childhood cancer diagnosis
model.dauno.ccage <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ Doxorubicin.Doseall.per100 + Daunorubicin.Doseall.per100
  + as.factor(Epirubicin.all) + as.factor(ChestRTfield.totaldose.com) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation)
  + as.factor(CED.cat) + Daunorubicin.Doseall.per100:as.factor(ccage.cat), data = IPD.nN.n0.2D, weights = weight)


# Additive interaction
#---------------------
# Interaction: Doxorubicin dose * Chest radiotherapy (yes/no)
model.doxo.chestRT.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Doxorubicin.Doseall.per100) * const(factor(ChestRT))
  + const(Daunorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ccage.cat)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)

# Interaction: Daunorubicin dose * Chest radiotherapy (yes/no)
model.dauno.chestRT.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Daunorubicin.Doseall.per100) * const(factor(ChestRT))
  + const(Doxorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ccage.cat)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)

# Interaction: Doxorubicin dose * Chest radiotherapy field/dose combination
model.doxo.chestRTfielddose.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Doxorubicin.Doseall.per100) * const(factor(ChestRTfield.totaldose.com))
  + const(Daunorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ccage.cat)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)

# Interaction: Daunorubicin dose * Chest radiotherapy field/dose combination
model.dauno.chestRTfielddose.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Daunorubicin.Doseall.per100) * const(factor(ChestRTfield.totaldose.com))
  + const(Doxorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ccage.cat)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)

# Interaction: Doxorubicin dose * Age at childhood cancer diagnosis
model.doxo.ccage.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Doxorubicin.Doseall.per100) * const(factor(ccage.cat))
  + const(Daunorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ChestRTfield.totaldose.com)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)

# Interaction: Daunorubicin dose * Age at childhood cancer diagnosis
model.dauno.ccage.add <- aalen(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ const(Daunorubicin.Doseall.per100) * const(factor(ccage.cat))
  + const(Doxorubicin.Doseall.per100) + const(factor(Epirubicin.all)) + const(factor(Pelvic5)) + const(factor(ChestRTfield.totaldose.com)) + strata(InstituteAbbreviation)
  + const(as.factor(CED.cat)), data = IPD.nN.nS.2D)


# Model: survivors with no Chest RT nor alklating agents use
#-----------------------------------------------------------
model.nochestRT.noalk <- coxph(Surv(time = entryage, time2 = exitage, event = SBC.DCIS.diag) ~ as.factor(Doxorubicin.cat.new) + as.factor(Daunorubicin.all)
    + as.factor(Epirubicin.all) + as.factor(Pelvic5) + as.factor(ccage.cat) + strata(InstituteAbbreviation), data = IPD.nN.n0.noc.noa.regression, weights = weight)


# Cumulative incidence
#---------------------
# by Chest radiotherapy (yes/no) * Doxorubicin dose categories
CI.Doxo.big.ChestRT <- survfit(Surv(time = entryage, time2 = exitage, event = status_ternary) ~ Doxorubicin.big.ChestRT, data = IPD.nN.n0, weights = weight, id = RecordId)
CI.Doxo.big.ChestRT

# by Doxorubicin dose categories in survivors with mantle field irradiation
CI.Doxo.mantle <- survfit(Surv(time = entryage, time2 = exitage, event = status_ternary) ~ Doxorubicin.cat.big, data = IPD.nN.n0.Mantle, weights = weight, id = RecordId)
CI.Doxo.mantle

# by Doxorubicin dose categories in survivors with mediastinal field irradiation
CI.Doxo.Mediastinal <- survfit(Surv(time = entryage, time2 = exitage, event = status_ternary) ~ Doxorubicin.cat.big, data = IPD.nN.n0.Mediastinal, weights = weight, id = RecordId)
CI.Doxo.Mediastinal

# by Doxorubicin dose categories in survivors with TBI or Whole lung field irradiation
CI.Doxo.TBI.Wholelung <- survfit(Surv(time = entryage, time2 = exitage, event = status_ternary) ~ Doxorubicin.cat.big, data = IPD.nN.n0.TBI.Wholelung, weights = weight, id = RecordId)
CI.Doxo.TBI.Wholelung

# by Doxorubicin dose categories in survivors with other chest fields irradiation
CI.Doxo.Otherfields <- survfit(Surv(time = entryage, time2 = exitage, event = status_ternary) ~ Doxorubicin.cat.big, data = IPD.nN.n0.Otherfields, weights = weight, id = RecordId)
CI.Doxo.Otherfields
