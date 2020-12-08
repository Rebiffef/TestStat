### EXO 5
# H0 : p=1/2 Vs H1 : p=p1
# binom.test

### EXO 6
# H0 : σ = 10 Vs H1 : σ = σ 1 < 10
X<-c(725, 722, 727, 718, 723, 731, 719, 724, 726, 726)

### EXO 7
# H0 : mu<=20 Vs H1 : mu>20
x<-c(20,23,23,23,22,20,23)
T_obs<-sqrt(length(x))*(mean(x)-20)/sd(x)
t_alpha<-qt(1-0.05,length(x)-1)
p_valeur<-1-pt(sqrt(length(x))*(mean(x)-20)/sd(x),length(x)-1)
t.test(x, mu = 20, alternative = "greater")

### EXO 9
# H0 : sigma²=10 Vs H1 : sigma²<>10
poids <- c(165.1,171.5,168.1,165.6,166.8,170.0,168.8, 171.1,168.8, 173.6,163.5,169.9,165.4,174.4,171.8,166.0,174.6,174.5,166.4,173.8)
n<-length(poids)
s_prime<-sd(poids)
T_obs<-(n-1)*s_prime^2/10
p_valeur<-2*min(pchisq(T_obs,n-1),1-pchisq(T_obs,n-1))
#install.packages("TRSbook")
#library(TRSbook)
### sigma 2
sigma2.test <- function (x, alternative = "two.sided", var0 = 1, conf.level = 0.95) {
  choices <- c("two.sided", "greater", "less")
  alt <- pmatch(alternative, choices)
  alternative <- choices[alt]
  if (length(conf.level) != 1 || is.na(conf.level) || conf.level <
      0 || conf.level > 1)
    stop("conf.level must be a number between 0 and 1")
  dname <- deparse(substitute(x))
  nx <- length(x)
  if (nx <= 2)
    stop("not enough x observations")
  gradiliberta <- nx - 1
  sx <- sd(x)
  estimate <- sx^2
  if (var0 > 0) s2obs <- (nx - 1) * sx^2/var0 else s2obs <- (nx - 1) * sx^2
  method <- c("One-sample Chi-squared test for given variance")
  names(estimate) <- c("var of x")
  if (var0 > 0) {
    if (alternative == "less") {
      pval <- pchisq(s2obs, df = nx - 1)
      cint <- c(0,(nx - 1) * sx^2/qchisq(p = 1 - conf.level,
                                         df = nx - 1))
    }
    else if (alternative == "greater") {
      pval <- 1 - pchisq(s2obs, df = nx - 1)
      cint <- c((nx - 1) * sx^2/qchisq(p = conf.level,
                                       df = nx - 1),Inf)
    }
    else {
      pval <- 2 * min(pchisq(s2obs, df = nx - 1), 1 - pchisq(s2obs,
                                                             df = nx - 1))
      cint <- c((nx - 1) * sx^2/qchisq(p = 1 - (1 - conf.level)/2,
                                       df = nx - 1), (nx - 1) * sx^2/qchisq(p = (1 - conf.level)/2,
                                                                            df = nx - 1))
    }
  } else {
    gradiliberta <- NA
    if (length(unique(x)) == 1) {
      pval <- 1
      cint <- c(0,0)
      alternative <- "greater"
    } else {
      pval <- 0
      cint <- c((nx - 1) * sx^2/qchisq(p = conf.level,
                                       df = nx - 1),Inf)
      alternative <- "greater"
    }
  }
  names(s2obs) <- "X-squared"
  names(gradiliberta) <- "df"
  names(var0) <- "variance"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = s2obs, parameter = gradiliberta,
               p.value = pval, conf.int = cint, estimate = estimate,
               null.value = var0, alternative = alternative, method = method,
               data.name = dname)
  attr(rval, "class") <- "htest"
  return(rval)
}
### sigma 2
sigma2.test(poids,var0=10)

### EXO 11
#H0 : p<=0.5 Vs H1:p>0.5 (point vue controleur)
#test exact 
citron<-c(6.71, 5.82, 6.16, 7.48, 6.43, 6.18, 6.79, 7.21, 7.74, 6.45, 6.82, 6.05, 6.51, 5.95, 6.13, 6.66, 6.79, 7.86, 5.93, 6.77, 7.32, 7.11, 6.95, 6.08, 6.65, 7.02, 6.63, 6.75, 6.78, 6.80, 6.44, 6.74, 7.22, 6.93, 6.18, 6.25, 6.42, 6.90, 7.03 , 7.24)
n<-length(citron)
p0<-0.5
T_obs<-length(citron[citron>6.5 & citron<7.3])
p_valeur<-1-pbinom(T_obs-1,n,p0)
binom.test(T_obs,n,p=0.5,alternative="greater") #=> on conserve H0
#test asymptotique
T_obs<-sqrt(40)*(22/40-0.5)/0.5
p_valeur<-1-pnorm(T_obs)
prop.test(22,40,alternative = "greater",correct=F) #=> on conserve H0

#TP 2 EXO 1
# hyp de modelisation : échantillons indépendants lois normales
# homoscedasticité
# H0 : sigma1=sigma2 Vs H1 : sigma1<>sigma2
M<-c(120, 107, 110, 116, 114, 111, 113, 117, 114, 112)
F<-c(110,111, 107, 108, 110, 105, 107, 106, 111, 111)
nM<-length(M)
nF<-length(F)
T_obs<-var(M)/var(F)
F_obs<-pf(T_obs,nM-1,nF-1)
p_val<-2*min(F_obs,1-F_obs) #-> on conserve H0
var.test(M,F)
# H0 : mu1=mu2 Vs H1 : mu1<>mu2
S<-((nM-1)*var(M)+(nF-1)*var(F))/(nF+nM-2)
t_obs<-(mean(M)-mean(F))/(sqrt(S*(1/nM+1/nF)))

F_obs<-pt(abs(t_obs),nM+nF-2)
p_valeur<-2*(1-F_obs)

t.test(M,F,var.equal = T,alternative="two.sided")

#TP2 EXO 2
# H0 : var(x)=var(y) Vs H1: var(x)<>var(y)
library(readr)
data <- read.csv("Intima_Media.csv",dec=",")

#data<_read.csv('C:/Users/fanny/Contacts/Desktop/test_stat/Intima_Media.csv')
#read.csv('C:\\Users\\fanny\\Contacts\\Desktop\\test_stat\\Intima_Media.csv')
x <- data[data$SPORT==1 & data$SEXE==2,]$mesure
y <- data[data$SPORT==0 & data$SEXE==2,]$mesure
n1<-length(x)
n2<-length(y)
T_obs<-var(x)/var(y)
F_obs<-pf(T_obs,n1-1,n2-1)
p_val<-2*min(F_obs,1-F_obs) #ccl : on rejette H0 -> var(x)<>var(y)

# TP2 exo 3
# H0 : mean(x)=mean(y)
t.test(x,y,var.equal = FALSE) # ccl : on conserve H0

#TP2 EXO 4
#H0 : pa=pb
na<-nb<-50
Fa<-27/50
Fb<-18/50
# on approche la loi de la différence des deux binomiales par une loi normale (n>50) car sinon on ne connait pas la loi
#H1 : pa<>pb
rho<-(na*Fa+nb*Fb)/(na+nb)
T_obs<-(Fa-Fb)/sqrt(rho*(1-rho)*(1/na+1/nb))
p_valeur<-2*(1-pnorm(abs(T_obs)))
prop.test(c(27,18),c(50,50),correct=F)
#H1 : pa>pb
prop.test(c(27,18),c(50,50),correct=F, alternative = "greater")
pval= 1-pnorm(Tobs)
# ccl : le test A est mieux que B