
setMethod(f="resid",signature=c(object="cold"), 

function (object, type = c( "pearson","response","null"),...)
{
    
	type <- match.arg(type)

	data<-object@subset.data

	data$y<-as.vector(t(object@y.matrix))  ### IMPORTANTE
	ynew<-data$y

        mu <- object@Fitted


        switch(type, pearson = , response = )


 	if(type=="pearson"||type=="null"){
   
	res = (ynew - mu)/sqrt(mu)}
   
	else
          {res = (ynew - mu)}

        names(res) <- c(1:length(res))
        round(res,6)
}
)

