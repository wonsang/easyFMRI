\name{meanROIs}
\alias{meanROIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The average time series of each ROI
}
\description{
Computing the average time series of each ROI in the given EPI sequence
}
\usage{
meanROIs(EPIfname, nROIs = 31, ROIdir = "", prefix = "roi_", 
     outfname = "meanROIs", writefile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EPIfname}{
the file name of a EPI sequence (The file format should be ANALYZE such as \code{hdr} and \code{img}. For example, if the original file is "\code{foo.img}", the input value should be "\code{foo}".)
}
  \item{nROIs}{
the number of ROIs. Default is \code{31}.
}
  \item{ROIdir}{
the directory of ROI definition files. These files should contain the individual zipped NIFTI files (nii.gz) corresponding to ROIs. Default is current working directory.
}
  \item{prefix}{
the common prefix of ROI definition files. Default is \code{roi_}. Then, the files should be titled as \code{roi_1.nii.gz}, \code{roi_2.nii.gz}, and so forth.
}
  \item{outfname}{
the outout text file name. Default is \code{meanROIs}.
}
  \item{writefile}{
Enable to write the output text file. Default is \code{TRUE}.
}
}
\details{
This function computes the average time series of each ROI in the given EPI sequence. The ROI definition files should be given a priori.
}
\value{
A matrix of averaged time series for ROIs. Row indicates a time point, and column does an ROI.
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
\code{\link{tsROI}}
}
\examples{
mts<-meanROIs("foo")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ averaged ROI }
\keyword{ FMRI }% __ONLY ONE__ keyword per line
