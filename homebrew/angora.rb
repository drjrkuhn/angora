# Documentation: https://github.com/Homebrew/brew/blob/master/docs/Formula-Cookbook.md
#                http://www.rubydoc.info/github/Homebrew/brew/master/Formula
# PLEASE REMOVE ALL GENERATED COMMENTS BEFORE SUBMITTING YOUR PULL REQUEST!

class Angora < Formula
  desc ""
  homepage ""
  url "http://www.angorafdtd.org/angora/angora-0.20.0.tar.gz"
  version "0.20.0"
  sha256 "67a8473e17c6bb278b6610177ee863bb4684c5de55b2d683669184cc81d8cbcf"

  depends_on :mpi => [:cc, :optional]

  depends_on "autoconf" => :build
  depends_on "automake" => :build
  depends_on "libtool" => :build
  depends_on "gettext" => :build

  mpi_args = Array.new
  mpi_args << "with-mpi" if build.with? "mpi"
  
  depends_on "hdf5" => mpi_args
  depends_on "h5utils" => :recommended
  depends_on "blitz"
  depends_on "libconfig"
  depends_on "argp-standalone"
  depends_on "boost" >> mpi_args
  

  def install
    conf_args = [
      "--disable-dependency-tracking",
      "--disable-silent-rules",
      "--prefix=#{prefix}"
    ]
    
    # export CXXFLAGS="-std=c++98 -Wno-parentheses"
    # export CPLUS_INCLUDE_PATH="/usr/local/include"
    # export LIBRARY_PATH="/usr/local/lib"
    # export LIBS="-largp -lhdf5 -lhdf5_cpp -lconfig++"
    # export HDF5_CXX="c++" HDF5_CLINKER="c++"

    ENV.append "CXXLAGS", "-I#{HOMEBREW_PREFIX}/include -std=c++98"
    ENV.append "LDFLAGS", "-L#{HOMEBREW_PREFIX}/lib -largp -lhdf5"

    # system "cmake", ".", *std_cmake_args
    system "make", "install" # if this fails, try separate make/make install steps
  end

  test do
    system "false"
  end
end
