setMethod("randeff",
          signature(object = "cold"),
          function (object) 
          {
            coef <- object@coefficients
            bi.est<- object@bi.estimate
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
            names(coef)<-cnames
            Z.aux<-object@random.matrix
            names.Z <- dimnames(Z.aux)[[2]]  #new for random
            X<-object@model.matrix
 
            # for one random effect
            if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) )
            {  n.cases<- length(bi.est)
          #     dimnames(bi.est) <- list(c(1:n.cases), "(Intercept)")
            id<-c(1:n.cases)
          #  tabela2<-data.frame(id,bi.est)
          #  names(tabela2)<-c(" ","(Intercept)")
          #  print(tabela2)
          #  tabela2<-data.frame(bi.est)
          #  names(tabela2)<-c("(Intercept)")
            tabela2<-as.data.frame(matrix(as.double(bi.est),nrow=n.cases,ncol=1))
            dimnames(tabela2) <- list(c(1:n.cases), "(Intercept)")
            return(tabela2)
            } 
            
            # for Two random effects
            else if (!all(is.na(match(cnames, "omega2"))) )
            { n.cases<- length(bi.est)/2
          #    dimnames(bi.est) <- list(c(1:n.cases), c("(Intercept)", names.Z[2]))
          #    print(bi.est)
            id<-c(1:n.cases)
          #  tabela2<-data.frame(id,bi.est[,1],bi.est[,2])
          #  names(tabela2)<-c(" ","(Intercept)", names.Z[2])
          #  print(tabela2)
            tabela2<-data.frame(bi.est[,1],bi.est[,2])
            names(tabela2)<-c("(Intercept)", names.Z[2])
            return(tabela2)
            }    
            
          else 
              warning("\nOnly to Random effects model") 
            
        })