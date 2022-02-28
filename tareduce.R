#useage:  At BASH prompt, type
#Rscript tareduce.R signalfilename backgroundfilename
#filenames should include the time and preceeding information, but not brightness.txt, which is added automatically
#tab autocompletion automatically helps you input the right file name info

#Laszlo Frazer, Borguet Group, 2015
#R version 3.1.1 (2014-07-10) -- "Sock it to Me"
#best to use this in a POSIX environment
#assumes gnuplot and sed are installed, and BASH is the shell
#assumes okular, a linux EPS viewer, is installed.  This is less important.



require('fBasics'); #for colSds
library(Rcpp); #for performance
library(foreach); #multicore operation
library(doMC);
#use all the cores
registerDoMC(detectCores());


#high performance functions
#remove zero and numbers below zero
cppFunction('List removezero(List inp) {
	    int length=inp.size();
	    for(int i=0;i<length;i++){
		    if(0>=inp[i]){
		    	inp[i]=NA_REAL;
		    }
	    }
	    return inp;
}');
cppFunction('List divide(List num,List den) {
	    List ans=num;
	    int length=num.size();
	    for(int i=0;i<length;i++){
		    ans[i]=(as<NumericVector>(num[i]))/(as<NumericVector>(den[i]));
	    }
	    return ans;
}');


#layout of data
pixels<-1600;
nondatarows<-5;

#get files to reduce from command line argument
#signal first, then background
file<-commandArgs(TRUE);
#file<-"sample/test2015-03-24-13-27-37"


#quit if wrong number of command line arguments
stopifnot(2==length(file));
cat("Signal File: ",file[1],"\nBackground File: ",file[2],"\n");

data<-read.table(paste(file[1],"brightness.txt",sep=""));
backdata<-read.table(paste(file[2],"brightness.txt",sep=""));

#read wavelength scale backwards since apparently CCD data is backwards
wavelength<-rev(read.table(paste(file[1],"wavelength.txt",sep="")));
backwavelength<-rev(read.table(paste(file[2],"wavelength.txt",sep="")));

#sanity check; if these are not the same we are using the wrong files
stopifnot(backwavelength==wavelength);

#subtract background, using mean value of background file
#assumes all pixels and times have the same background
print("Subtracting Background");
#there are 5 columns of metadata (nondata rows) we do not need to subtract
data[,(1:pixels)+nondatarows]<-data[,(1:pixels)+nondatarows]-mean(as.matrix(backdata[,(1:pixels)+nondatarows]));


print("Computing OD");

colnames(data)<-c("positionnumber","readnumber","positionmm","positiont","phase");


#ensure even number of rows 
stopifnot(0==nrow(data)%%2);
#ensure there are at least four rows
stopifnot(3<nrow(data));

#initialize dataframes for speed
positionnumber<-0;
readnumber<-0;
OD<-data.frame(matrix(ncol=pixels+1,nrow=data[nrow(data),1]+1));
ODsd<-data.frame(matrix(ncol=pixels+1,nrow=data[nrow(data),1]+1));
trans<-data.frame(matrix(ncol = pixels, nrow = 
	(nrow(subset(subset(data,positionnumber==0),readnumber==0))%/%2-1)
));

for(i in 0:data[nrow(data),1]){
	cat("Position: ",i,"\n");
	
	#select data for this position
	pdata<-subset(data,positionnumber==i);

	#initialize ODs for this position
	#TODO:  unroll this
	ptrans<-{};
	
	#loop through all reads in parallel
	ptrans<-foreach(j = 0:data[nrow(data),2],.combine='rbind')%dopar%{

		#select data for this read
		rdata<-subset(pdata,readnumber==j);

		#ensure even number of rows 
		stopifnot(0==nrow(rdata)%%2);

		#skip first pair of rows in each read because first row is contaminated by lack of keep cleans
		for(k in 2:(nrow(rdata)%/%2)){

			#determine phase, proportion transmitted
			#may produce Inf or NaN for excessively weak probe
			if(1==rdata[1,]$phase){
				#pulse 2k-1 has pump on
				
				#replace Inf/NaN with NA
				#by removing zero and negtaive values from the denomenator
				#there are 5 columns of metadata we do not need to include
				denomenator<-removezero(rdata[2*k,(1:pixels)+nondatarows]);

				#compute transmission
				trans[k-1,]<-divide(rdata[2*k-1,(1:pixels)+nondatarows],denomenator);

			} else if (0==rdata[1,]$phase){
				#pulse 2k-1 has pump off
				#so use reciprocal
				
				#replace Inf/NaN with NA
				#by removing zero and negtaive values from the denomenator
				#there are 5 columns of metadata we do not need to include
				denomenator<-removezero(rdata[2*k-1,(1:pixels)+nondatarows]);

				#compute transmission
				trans[k-1,]<-divide(rdata[2*k,(1:pixels)+nondatarows],denomenator);
			} else {
				for(l in (1:pixels)+nondatarows){
					trans[k-1,l-5]<-NA;
				}
			}
		}
		trans;
	}


	#average OD for this stage position
	#compute log for OD after averaging so that negative outliers get averaged away
	#include delay time in first column
	#this will produce and propagate -Inf and NaN for excessively weak probe
	OD[i+1,]<-cbind(pdata[1,]$positiont,t(-log10(colMeans(ptrans,na.rm=TRUE))));
	
	ODsd[i+1,]<-cbind(pdata[1,]$positiont,t(
						#uncertainty in log10 x=  abs(dx/(x*ln(10)))
						abs(colSds(ptrans,na.rm=TRUE)/(colMeans(ptrans,na.rm=TRUE)*log(10)))
						));
}

#formatting first row for gnuplot
firstrow<-cbind(pixels+1,wavelength);

#strip conflicting names
names(OD)<-rep("name",ncol(OD));
names(ODsd)<-rep("name",ncol(ODsd));
names(firstrow)<-rep("name",ncol(OD));

#write data to files
write.table(rbind(firstrow,OD),paste(file[1],"reduced.tsv",sep=""),row.names=FALSE,col.names=FALSE);
write.table(rbind(firstrow,ODsd),paste(file[1],"stdev.tsv",sep=""),row.names=FALSE,col.names=FALSE);
print("Analysis Done");

#copy the gnuplot template, inserting correct file names
system(paste("sed 's/PLACEHOLDER/",file[1],"/' plottemplate.gp > ",file[1],".gp",sep=""));

#run gnuplot 
system(paste("gnuplot -e \"load '",file[1],".gp'\"",sep=""));

#display graph
system(paste("okular ",file[1],".eps &",sep=""));
system(paste("okular ",file[1],"line.eps &",sep=""));
print("Plot Done");

print(warnings());
