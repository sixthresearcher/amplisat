##############################################################################################
#
# Dockerfile for AmpliSAT (Amplicon Sequencing Analysis Tools)
# 
# You can use the tools online at:
#   http://evobiolab.biol.amu.edu.pl/amplisat
#
# Author: Alvaro Sebastian
#   http://www.sixthresearcher.com/
#
##############################################################################################

# Based on
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Alvaro Sebastian <sixthresearcher@gmail.com>

# Update repos and install basic software that could be required later
RUN apt-get update \
 && apt-get install -y --no-install-recommends software-properties-common build-essential cpanminus expat libexpat1-dev libipc-run-perl libssl-dev \
 && apt-get install -y --no-install-recommends git vim less wget curl unzip file emboss

# Install required Perl modules
RUN cpanm Bio::Seq \
 && cpanm IO Sort::Naturally Archive::Zip Excel::Writer::XLSX Spreadsheet::XLSX \
 && cpanm Text::Iconv Net::SFTP::Foreign Statistics::Descriptive \
 && cpanm File::FindLib Inline::C threads
 
# Clone the AmpliSAT repo and tests if it works
RUN cd /opt \
 && git clone https://github.com/sixthresearcher/amplisat.git \
 && echo $(/opt/amplisat/ampliSAS.pl -h)
ENV PATH="/opt/amplisat:${PATH}"

# Set up a working dir
WORKDIR /workdir




