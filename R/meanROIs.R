meanROIs <-
function(EPIfname, nROIs=31, ROIdir="", prefix="roi_", outfname="meanROIs", writefile=TRUE) {
	
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
	mtsROIs<-mat.or.vec(nROIs, nTPs)
	for(i in seq(1,nROIs)) {
		ROIfname<-paste(c(ROIdir,prefix,i,".nii.gz"),sep="",collapse="")
		roi<-readNIfTI(ROIfname)@.Data
		roipos<-roi>0
		for(t in seq(1,nTPs)){
			#epiTP<-epi[,,,t]
			epiTP<-f.read.nifti.tpt(EPIimg,t)
			epiTPROI<-epiTP[roipos]
			mtsROIs[i,t]<-mean(epiTPROI)
		}
	}
	
	mtsROIs
	
	# Write a text file
	if(writefile)
		write.table(mtsROIs, file = paste(c(outfname,".txt"),sep="",collapse=""), 
				sep = ";", col.names = NA)
}

