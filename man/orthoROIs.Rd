\name{orthoROIs}
\alias{orthoROIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonalization of ROIs
}
\description{
Orthogonalizes ROIs in the EPI sequence. In other words, it eliminates all overlapped parts among ROIs, and then creates new ROIs which are spatially independent each other.
}
\usage{
orthoROIs(nROIs = 31, ThRes = 5, ROIdir = "", 
      targetdir = "newROIs", mode = "minimal", prefix = "roi_", 
      newprefix = "roi_", writefile = TRUE, gzipped = TRUE, 
      verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nROIs}{
the number of ROIs
}
  \item{ThRes}{
the resolution of threshold. It is applied only when you choose the mode as thresholding. For example, if \code{ThRes=5}, the precision becomes \code{1e-5}. Default is \code{5}. 
}
  \item{ROIdir}{
the directory of ROI definition files. These files should contain the individual zipped NIFTI files (nii.gz) corresponding to ROIs. Default is current working directory.
}
  \item{targetdir}{
the directory of output ROI definition files. Default is \code{newROIs}.
}
  \item{mode}{
the mode of orthogonalization. \code{thres} denotes thresholding, and \code{minimal} denotes minimal overlapping elimination. Default is \code{minimal}.
}
  \item{prefix}{
the common prefix of ROI definition files. Default is \code{roi_}. Then, the files should be titled as \code{roi_1.nii.gz}, \code{roi_2.nii.gz}, and so forth.
}
  \item{newprefix}{
the common prefix of output ROI definition files. Default is \code{roi_}.
}
  \item{writefile}{
Enable to write the output text file. Default is \code{TRUE}.
}
  \item{gzipped}{
Enable out NIFTI file to be zipped. Default is \code{TRUE}.
}
  \item{verbose}{
Allows message notification duringexecution of the function. 
}
}
\details{
ROIs is normally defined from the anatomical high-resolution image of the brain. To define the ROIs in the EPI sequence, the Brain Voyager can be exploited, however it can cause seriosu distortion such as spatial and temporal blurring. To avoid temporal blurring, we can match each ROI into the EPI space through rotation, translation, and resampling. In this case, overlaps among ROIs can be generated. This function resizes all ROIs by eliminating all overlapped parts.
}
\value{
  \item{ROIData}{
A list of matrices of ROI data. Each matrix is three-dimensional, and indicates a segmented ROI region at the 3D volume of FMRI data. The number of matrices is the same as the number of ROIs.
}
  \item{thr}{
A vector of thresholds of ROIs. Notice that each ROI can have different threshold. The length of this vector is the same as the number of ROIs.
}
  \item{nOverlapROIpairs}{
An integer which indicates the number of overlapped ROI pairs.
}
  \item{nOverlapVoxels}{
An integer which indicates the number of overalpped voxels.
}
  \item{nNewOverlapROIpairs}{
An integer which indicates the number of overlapped ROI pairs after orthogonalization. In general, it is expected to be zero.
}
  \item{nNewOverlapVoxels}{
An integer which indicates the number of overalpped voxels after orthogonalization. In general, it is expected to be zero.
}
}
\references{
Wonsang You (2010) ROI Data Extraction from FMRI BOLD Signals of the Human Brain, Technical Reports of the Leibniz Institute for Neurobiology, TR 10016.
}
\author{
Wonsang You
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dat<-orthoROIs()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ROI orthogonalization }
\keyword{ ROI overlap }% __ONLY ONE__ keyword per line
