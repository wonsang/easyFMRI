orthoROIs <-
function(nROIs=31, ThRes=5, ROIdir="", targetdir="newROIs", mode="minimal", prefix="roi_", newprefix="roi_", writefile=TRUE, gzipped=TRUE, verbose=TRUE) {
	
	require(oro.nifti)
	
	# Extract ROI data
	if(length(ROIdir)>1)
		ROIdir<-paste(c(ROIdir,"/"), collapse="")
	
	ROIdata<-c()
	for(i in seq(1,nROIs)) {
		roi<-readNIfTI(paste(c(ROIdir,prefix,i,".nii.gz"),sep="",collapse=""))@.Data
		ROIdata<-c(ROIdata,list(roi))
	}
	
	# Evaluate overlaps
	olMat<-mat.or.vec(nROIs, nROIs)
	thr<-mat.or.vec(nROIs,1)
	for(i in seq(1,nROIs-1)) {
		x1<-ROIdata[[i]]
		olMat[i,i]<-0
		for(j in seq(i+1,nROIs)) {			
			x2<-ROIdata[[j]]
			olMat[i,j]<-length(x1[x1&x2])
		}
	}
	
	# The number of overlapped ROIs
	nOverlapROIpairs<-length(olMat[olMat>0])
	nTotalROIpairs<-nROIs*(nROIs-1)/2
	
	# The number of overlapped Voxels
	nOverlapVoxels<-sum(olMat[olMat>0])

	if(verbose){
		msg <- paste(c("The number of overlapped ROI pairs", nOverlapROIpairs,"/",nTotalROIpairs),collapse = " ")
		print(msg)
		msg <- paste(c(" The number of overlapped Voxels= ", nOverlapVoxels),collapse = " ")
		print(msg)
	}
	
	# Resize ROIs to eliminate overlaps
	olMat<-mat.or.vec(nROIs, nROIs)
	thr<-mat.or.vec(nROIs,1)
	for(i in seq(1,nROIs-1)) {
		x1<-ROIdata[[i]]
		olMat[i,i]<-0
		for(j in seq(i+1,nROIs)) {			
			x2<-ROIdata[[j]]
			olMat[i,j]<-length(x1[x1&x2])
			
			if(olMat[i,j]>0){
				if(mode=="minimal"){
					x1<-replace(x1,(x1>0)&(x2>0),0)
					ROIdata[[j]]<-replace(x2,(x1>0)&(x2>0),0)
				}else{
					nOverlapVoxelsInROI<-olMat[i,j]
					
					ux1<-x1
					ux2<-x2
					thr1<-thr[i]
					thr2<-thr[j]
					for(r in seq(1,trunc(ThRes))){
						ThrRes<-1/10^r
						while(nOverlapVoxelsInROI>0){
							tmp_thr1<-thr1+ThrRes
							tmp_thr2<-thr2+ThrRes
							tmp_ux1<-replace(ux1,ux1<thr1,0)
							tmp_ux2<-replace(ux2,ux2<thr2,0)
							nOverlapVoxelsInROI<-length(tmp_ux1[tmp_ux1&tmp_ux2])
							
							if(nOverlapVoxelsInROI>=0){
								thr1<-tmp_thr1
								thr2<-tmp_thr2
								ux1<-tmp_ux1
								ux2<-tmp_ux2
							}
						}
					}
					x1<-ux1
					ROIdata[[j]]<-ux2
					thr[i]<-thr1
					thr[j]<-thr2	
				}			
			}
		}
		ROIdata[[i]]<-x1
	}
	
	# Verify overlaps
	olMat<-mat.or.vec(nROIs, nROIs)
	for(i in seq(1,nROIs-1)) {
		x1<-ROIdata[[i]]
		olMat[i,i]<-0
		for(j in seq(i+1,nROIs)) {			
			x2<-ROIdata[[j]]
			olMat[i,j]<-length(x1[x1&x2])
		}
	}
	
	nNewOverlapROIpairs<-length(olMat[olMat>0]) # The number of overlapped ROIs
	nNewOverlapVoxels<-sum(olMat[olMat>0]) # The number of overlapped Voxels
	
	if(verbose){
		msg <- paste(c("The number of new overlapped ROI pairs", nNewOverlapROIpairs,"/",nTotalROIpairs),collapse = " ")
		print(msg)
		msg <- paste(c(" The number of new overlapped Voxels= ", nNewOverlapVoxels),collapse = " ")
		print(msg)		
	}
	
	# Write NIFTI files
	if(writefile){
		dir.create(paste(c(ROIdir,targetdir),sep="",collapse=""))
		for(i in seq(1,nROIs)){
			roi.fname<-paste(c(ROIdir,targetdir,"/",newprefix,i),sep="",collapse="")
			roi.nifti<-nifti(ROIdata[[i]],datatype=64) # double type
			writeNIfTI(roi.nifti, roi.fname, gzipped=gzipped, verbose=FALSE)
		}
	}
	
	if(mode=="minimal"){
		dat<-drop(list(ROIData=ROIdata,nOverlapROIpairs,nOverlapVoxels,
						nNewOverlapROIpairs,nNewOverlapVoxels))
	}else{
		dat<-drop(list(ROIData=ROIdata,thr=thr,nOverlapROIpairs,nOverlapVoxels,
						nNewOverlapROIpairs,nNewOverlapVoxels))
	}
	invisible(dat)
}

