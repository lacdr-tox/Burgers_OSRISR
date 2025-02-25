# Model for compound dynamics scaling = 100 (NIT) or 2500 (DIC&KET)
Compound_M1 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    list(c(dS, dCompound))
  })
}

# Adapted model for compound dynamics for ketoconazole
Compound_M2 <- function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*(Compound**n)
    list(c(dS, dCompound))
  })
}

# Model for OSR and ISR induction
# Separate pathway induction, no cross-talk, full model
Pathways_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade * Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = buildCBase + Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = buildK1Base - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1) - degradK1 * KEAP1
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) - degradK1 * modified_K1
    dNRF2 = buildN2Base - Vmax_degN2*((KEAP1*NRF2)/(K_N2 + NRF2)) - degradN2* NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Separate pathway induction, no cross-talk, simplified model
Pathways_M2 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Separate pathway induction, NRF2 transcription via ATF4, simplified model
Pathways_M3 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Separate pathway induction, NRF2 transcription via ATF4 withut bound, full model
Pathways_M4 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade * Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2* NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Model M3, adapted for ketaconazole
# added sc_OSR & sc_ISR, different function for compound
PathwaysKET_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*(Compound**n)
    dATF4 = buildA4Base + sc_ISR * V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (sc_OSR * k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (sc_OSR * k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Model M3 adapted for diclofenac
# added sc_OSR & sc_ISR
PathwaysDIC_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + sc_ISR * V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (sc_OSR * k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (sc_OSR * k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1))
  })
}

# Basic model for GSH 
GSH_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 =  - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2* NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    dGSH = buildGBase + Vmax_G * ((NRF2**nG)/(K_G**nG + NRF2**nG)) - degradG * GSH
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1, dGSH))
  })
}

# Model for GSH with ATF4-induced GSH degradation
GSH_M2 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 =  - (k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2* NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    dGSH = buildGBase + Vmax_G * ((NRF2**nG)/(K_G**nG + NRF2**nG)) - degradG * (degGA4 * ATF4 + 1) * GSH
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1, dGSH))
  })
}

# GSH Model 2 of nitrofurantoin adapted for diclofenac
GSHDIC_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*Compound
    dATF4 = buildA4Base + sc_ISR * V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (sc_OSR * k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (sc_OSR * k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    dGSH = buildGBase + Vmax_G * ((NRF2**nG)/(K_G**nG + NRF2**nG)) - degradG * (degGA4 * ATF4 + 1) * GSH
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1, dGSH))
  })
}

# GSH Model 2 of nitrofurantoin adapted for ketoconazole
GSHKET_M1 = function(t, inistate, parameters) {
  with(as.list(c(inistate, parameters)), {
    dS = -tau1*S
    dCompound = S - degrade*(Compound**n)
    dATF4 = buildA4Base + sc_ISR * V_A4S * Compound - degradA4 * ATF4
    dCHOP = Vmax_C * ((ATF4**nC)/(K_C**nC + ATF4**nC)) - degradC * CHOP
    dKEAP1 = - (sc_OSR * k_K1_mod*Compound*KEAP1) + (k_K1_unmod*modified_K1)
    dmodified_K1 = (sc_OSR * k_K1_mod*Compound*KEAP1) - (k_K1_unmod*modified_K1) 
    dNRF2 = buildN2Base + V_N2A4*ATF4 - V_degN2*(KEAP1*NRF2) - degradN2 * NRF2
    dSRXN1 = buildS1Base + Vmax_S1 * ((NRF2**nS1)/(K_S1**nS1 + NRF2**nS1)) - degradS1 * SRXN1
    dGSH = buildGBase + Vmax_G * ((NRF2**nG)/(K_G**nG + NRF2**nG)) - degradG * (degGA4 * ATF4 + 1) * GSH
    list(c(dS, dCompound, dATF4, dCHOP, dKEAP1, dmodified_K1, dNRF2, dSRXN1, dGSH))
  })
}