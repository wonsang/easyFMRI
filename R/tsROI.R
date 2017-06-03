tsROI <-
function(ROIidx, EPIfname, ROIdir="", prefix="roi_", outfname="tsROI", writefile=TRUE) {
	
	require(oro.nifti)
	require(AnalyzeFMRI)
	
	# Extract ROI data
	if(length(ROIdir)>1)
		ROIdir<-paste(c(ROIdir,"/"), collapse="")
	
	# Read the EPI sequence
	#dimEPI<-dim(EPIdata)
	#nTPs<-dimEPI[4]
	EPIimg<-paste(c(EPIfname,".img"),sep="",collapse="")
	EPIhdr<-paste(c(EPIfname,".hdr"),sep="",collapse="")
	nTPs<-f.read.nifti.header(EPIhdr)$dim[5]
	
	# Compute mean time series of each ROI
	ROIfname<-paste(c(ROIdir,prefix,ROIidx,".nii.gz"),sep="",collapse="")
	roi<-readNIfTI(ROIfname)@.Data
	roipos<-roi>0
	ts<-c()
	for(t in seq(1,nTPs)){		
		epiTP<-f.read.nifti.tpt(EPIimg,t)
		#epiTP<-EPIdata[,,,t]
		epiTPROI<-epiTP[roipos]
		ts<-rbind(ts,epiTPROI)
	}

	ts
	
	# Write a text file
	if(writefile)
		write.table(ts, file = paste(c(outfname,".txt"),sep="",collapse=""), 
				sep = ";", col.names = NA)
}

