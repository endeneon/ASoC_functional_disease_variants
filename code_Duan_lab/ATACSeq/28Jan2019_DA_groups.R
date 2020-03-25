# Siwei 28 Jan 2019
# peak=summit+/-250bp
# Identify DA peaks for each group (NPC_iPS, NPC_Glut, NPC_GA, NPC_DN)
# Use Xin's original summit set

library(DiffBind)
library(gplots)
library(factoextra)
library(grDevices)
library(colorRamps)
library(ggplot2)

### NPC vs iPS

NPC_iPS_contrast <- dba.contrast(DiffBind_read_counted_summit, # use iPS as baseline since NPC is differentiated from iPS
                                 group1 = DiffBind_read_counted_summit$masks$neural_progenitor_cell, 
                                 group2 = DiffBind_read_counted_summit$masks$induced_pluripotent, 
                                 name1 = "NPC", name2 = "iPS")
NPC_iPS_contrast <- dba.analyze(NPC_iPS_contrast, method = DBA_EDGER, bTagwise = T, 
                                bReduceObjects = F, bParallel = T)
NPC_iPS_contrast$config$AnalysisMethod <- DBA_EDGER
NPC_iPS_DA_peakset <- dba.report(NPC_iPS_contrast, method = DBA_EDGER, 
                                 th = 0.001, bUsePval = F)
plot(NPC_iPS_contrast, method = DBA_EDGER_GLM, contrast = 1)

### Glut vs NPC

Glut_NPC_contrast <- dba.contrast(DiffBind_read_counted_summit, # use iPS as baseline since NPC is differentiated from iPS
                                 group1 = DiffBind_read_counted_summit$masks$cortical, 
                                 group2 = DiffBind_read_counted_summit$masks$neural_progenitor_cell, 
                                 name1 = "Glut", name2 = "NPC")
Glut_NPC_contrast <- dba.analyze(Glut_NPC_contrast, method = DBA_EDGER, bTagwise = T, 
                                bReduceObjects = T, bParallel = T)
# NPC_iPS_contrast$config$AnalysisMethod <- DBA_EDGER
Glut_NPC_DA_peakset <- dba.report(Glut_NPC_contrast, method = DBA_EDGER, 
                                 th = 0.001, bUsePval = F)
plot(NPC_iPS_contrast, method = DBA_EDGER, contrast = 1)

#### GA vs NPC

GA_NPC_contrast <- dba.contrast(DiffBind_read_counted_summit, # use iPS as baseline since NPC is differentiated from iPS
                                  group1 = DiffBind_read_counted_summit$masks$GABAergic, 
                                  group2 = DiffBind_read_counted_summit$masks$neural_progenitor_cell, 
                                  name1 = "GABAergic", name2 = "NPC")
GA_NPC_contrast <- dba.analyze(GA_NPC_contrast, method = DBA_EDGER, bTagwise = T, 
                                 bReduceObjects = T, bParallel = T)

GA_NPC_DA_peakset <- dba.report(GA_NPC_contrast, method = DBA_EDGER, 
                                  th = 0.001, bUsePval = F)
#### DN vs NPC

DN_NPC_contrast <- dba.contrast(DiffBind_read_counted_summit, # use iPS as baseline since NPC is differentiated from iPS
                                group1 = DiffBind_read_counted_summit$masks$dopaminergic, 
                                group2 = DiffBind_read_counted_summit$masks$neural_progenitor_cell, 
                                name1 = "Dopaminergic", name2 = "NPC")
DN_NPC_contrast <- dba.analyze(DN_NPC_contrast, method = DBA_EDGER, bTagwise = T, 
                               bReduceObjects = F, bParallel = T)

DN_NPC_DA_peakset <- dba.report(DN_NPC_contrast, method = DBA_EDGER, 
                                th = 0.001, bUsePval = F)

############################

