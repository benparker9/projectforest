#install packages
install.packages('dplyr')
install.packages('tidyr')
install.packages('ggplot2')
#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
#define path to yield table, sort, clean data
dats <- read.csv("c:Users/xxx/YIELD.csv", header = TRUE)
datas <- dats %>%
  dplyr::group_by(Latin_Name)%>%
  dplyr::select(Latin_Name, Age, DBH)%>%
  arrange(Latin_Name)
print(datas, n=500)
#create output for loop
save <- c()
qul <- c()
#model for loop
for (i in unique(datas$Latin_Name)){
  filt <- datas %>%
    filter(Latin_Name == i)%>% dplyr::select(Age, DBH)
  #define cycle of species names
  is <- i
  
  #2nd order poly fitting
  poly <- lm(filt$DBH ~ poly(filt$Age, 2, raw = TRUE), data =filt)
  GoF2 <- summary(poly)$r.squared
  A <- coef(poly)[3]
  B <- coef(poly)[2]
  C <- coef(poly)[1]
  poly2 <- function(x){(A*x^2)+(B*x)+C}
  equation2 <- sprintf("(%s*Age^2) + (%s*Age) + %s", A, B, C)
  pred2 <- poly2(50)
  sum <- data.frame(x=1:(2*max(filt$Age)), y=poly2(1:(2*max(filt$Age))))
  ex2 <- poly2(1:50)
  #3rd order polynomial fit
  poly3 <- lm(DBH ~ poly(Age, 3, raw = TRUE), data =filt)
  GoF3 <- summary(poly3)$r.squared
  A3 <- coef(poly3)[4]
  B3 <- coef(poly3)[3]
  C3 <- coef(poly3)[2]
  D3 <- coef(poly3)[1]
  polyf3 <- function(x){(A3*x^3) + (B3*x^2) + (C3 * x) + D3}
  equation3 <- sprintf("(%s*Age^3) + (%s*Age^2) + (%s * Age) + %s", A3, B3, C3, D3)
  pred3 <- polyf3(50)
  curvepoly3 <- data_frame(Year = 1: (2*max(filt$Age)), DBH = polyf3(1:(2*max(filt$Age))))
  ex3 <- polyf3(1:50)
  #power function
  pwr <- nls(DBH ~ A * Age ^B, 
             data = filt,
             start = list(A=10, B=.5),
             control = list(maxiter = 1000),
             trace = TRUE,
             alg = "port")
  GoFpwr <- cor(filt$DBH,predict(pwr))
  Apwr <- coef(pwr)[1]
  Bpwr <- coef(pwr)[2]
  pwrf <- function(x){Apwr * x ^Bpwr}
  equationpwr <- sprintf("%s * Age ^ %s", Apwr, Bpwr)
  predpwr <- pwrf(50)
  curvepwr <- data.frame(Year = 1: (2*max(filt$Age)), DBH = pwrf(1:(2*max(filt$Age))))
  expwr <- pwrf(1:50)
  #log function
  log <- nls(DBH ~ A + B*log(Age), 
             data = filt,
             start = list(A=10, B=.5),
             control = list(maxiter = 1000),
             trace = TRUE,
             alg = "port")
  GoFlog <- cor(filt$DBH,predict(log))
  Alog <- coef(log)[1]
  Blog <- coef(log)[2]
  logf <- function(x){Alog + Blog * log(x)}
  predlog <- logf(50)
  equationlog <- sprintf("%s + %s * log(Age)", Alog, Blog)
  curvelog <- data.frame(Year = 1: (2*max(filt$Age)), DBH = logf(1:(2*max(filt$Age))))
  exlog <- logf(1:50)
  
  #cr fit
  tryCatch ({
    cr <- nls(DBH ~ A *(1 - exp(-K * Age))^P, 
              data = filt,
              start = list(A=max(filt$DBH), K=0.03, P=1),
              control = list(maxiter = 1000),
              trace = TRUE,
              alg = "port")
    GoFcr <- cor(filt$DBH,predict(cr))
    GoFcr
    Acr <- coef(cr)[1]
    Bcr <- coef(cr)[2]
    Ccr <- coef(cr)[3]
    crf <- function(x){Acr * (1-exp(-Bcr * x)) ^ Ccr}
    predcr <- crf(50) 
    equationcr <- sprintf("%s * (1-exp(-%s * Age)) ^ %s", Acr, Bcr, Ccr)
    curvecr <- data.frame(Year = 1:max(filt$Age), DBH = crf(1:max(filt$Age)))
    excr <- crf(1:50)
  },
  error = function(msg){
    return(NA)})
  
#spreadsheet to update
  mas <- data.frame(Species = i,
                    Poly2DBH = ex2,
                    Poly3DBH = ex3,
                    PwrDBH = expwr,
                    LogDBH = exlog,
                    CRDBH = excr)
  save <- append(save, mas)
  
#create comparison tool
  ql <- data.frame(Species = i,
                   DBH = c(pred2, pred3, predpwr, predlog, predcr), 
                   R2 = c(GoF2, GoF3, GoFpwr, GoFlog, GoFcr), 
                   Fit = c(equation2, equation3, equationpwr, equationlog, equationcr))
  
  qul<- append(qul, ql)
  
  #vizualise
  viz <- ggplot(data = filt, aes(x=Age, y=DBH))+
    geom_point()+
    geom_line(data=sum, aes(x=x, y=y, colour="2nd Order Polynomial"))+
    geom_line(data=curvepoly3, aes(x=Year, y=DBH, color="3nd Order Polynomial Prediction"))+
    geom_line(data=curvepwr, aes(x=Year, y=DBH, color="Power Function"))+
    geom_line(data=curvelog, aes(x=Year, y=DBH, color="Logarithmic Function"))+
    geom_line(data=curvecr, aes(x=Year, y=DBH, color="Chapman Richard's Prediction"))+
    labs(title=sprintf("Life Cycle of %s", is),
         x="Age (Years)", y="DBH (cm)", 
         subtitle="Data Gathered by RRG")+
    xlim(c(0,(2*max(filt$Age))))+
    ylim(c(0,(2*max(filt$DBH))))+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))+
    theme_light()
  #print viz and quick look
  print(viz)
  print(ql)
  #export
  forms <- sprintf("C:/Users/xxx/lifespan%s.jpg", is)
  ggsave(viz, file=forms)
}
#save and export tables
form <- "C:/Users/xxx/MASTERPredict.csv"
write.csv(save, file = form)
qulf <- "C:/Users/xxx/MASTERStatResults.csv"
write.csv(qul, file=qulf)