///////////////////////////////////////////////////////////////////////////////
// Convert several pgm image into a single 3dpgm file.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <deque>


#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"



using namespace std;
using namespace ImaGene;

static Arguments args;



struct ImagePGM{
  uint width;
  uint height;
  uint maxValue;
  unsigned char  *tabImage;
};


struct ImagePGM3D{
  uint width;
  uint height;
  uint nbSlice;
  uint maxValue;
  unsigned char  *tabImage;
};



ImagePGM3D importPGM3D(  istream &in);
ImagePGM importPGM(const char * name);
ImagePGM3D importPGM3D(const  char * name);
void exportPGM3D(ostream &out, ImagePGM3D image);


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addDigitalArgs( args, 2, false, false );
  args.addOption( "-init", "-init <image.pgm>: init an initial pgm3d file based on the first pgm image.", " " );
  args.addOption( "-addImage", "-addImage <image.pgm>: add to the std input the image as a new slice (and upate pgm3d header )", " " );



  

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) )     {
      cerr << args.usage( "pgm2pgm3D", 
			  "Simple tools to convert pgm files into a single pgm3d file",
			  "-init -addImage" )
	   << endl;
      return 1;
    }




  if(args.check("-init")){
    string nomInitialImage = args.getOption("-init")->getValue(0);
    const char * name = nomInitialImage.c_str();
    ImagePGM img = importPGM(name);  
    ImagePGM3D img3d;
    img3d.width = img.width;
    img3d.height = img.height;
    img3d.nbSlice = 1;
    img3d.maxValue = img.maxValue;
    img3d.tabImage = img.tabImage;
    exportPGM3D(cout, img3d);
    delete img.tabImage;
    return 0;
  }
  

  if(args.check("-addImage")){
    ImagePGM3D img3d = importPGM3D(cin); 
    const char * name = (args.getOption("-addImage")->getValue(0)).c_str();
    ImagePGM img = importPGM(name);
    uint tailleIni =img3d.width*img3d.height*(img3d.nbSlice); 
    img3d.nbSlice++;
    unsigned char * newTab = new unsigned char[tailleIni+img3d.width*img3d.height];
    // adding the new slice into the new tab
    for(int i=0; i< tailleIni ; i++){
      newTab[i] = img3d.tabImage[i];      
    }
    
    for(int i=0; i<img.width*img.height; i++){
      newTab[i+tailleIni]= img.tabImage[i];      
    }
    delete img3d.tabImage;
    img3d.tabImage= newTab;
    exportPGM3D(cout, img3d);
    delete img.tabImage;
    delete newTab;
  }  
  
  
  return 0;
}






ImagePGM 
importPGM(const  char * name){
  ifstream in;
  in.open(name, ifstream::in);
  
  ImagePGM imageImported;
  string str;
  getline( in, str );
  if ( ! in.good() ) return imageImported;
  if ( str != "P5" ) return imageImported;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return imageImported;
    }
  while ( str[ 0 ] == '#' || str=="" );
  istringstream str_in( str );
  str_in >> imageImported.width >> imageImported.height;
  in >> noskipws;
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  
  str2_in >> max_value;
  imageImported.maxValue=max_value;
  imageImported.tabImage = new unsigned char [imageImported.width*imageImported.height];  
  
  for(int i=0; i<imageImported.width*imageImported.height; i++){
    unsigned char c; 
    in >> c;
    imageImported.tabImage[i] = c;        
  }  
  return imageImported;
}







ImagePGM3D 
importPGM3D(const  char * name){
  ifstream in;
  in.open(name, ifstream::in);
  return importPGM3D(in);
}





ImagePGM3D
importPGM3D(istream &in){
  ImagePGM3D imageImported;
  string str;
  getline( in, str );
  if ( ! in.good() ) return imageImported;
  if ( str != "P3d" ) return imageImported;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return imageImported;
    }
  while ( str[ 0 ] == '#' || str == "" );
  istringstream str_in( str );
  str_in >> imageImported.width >> imageImported.height >> imageImported.nbSlice;
  in >> noskipws;
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  
  str2_in >> max_value;
  imageImported.maxValue=max_value;
  imageImported.tabImage = new unsigned char [imageImported.width*imageImported.height*imageImported.nbSlice];  
  
  for(int i=0; i<imageImported.width*imageImported.height*imageImported.nbSlice; i++){
    unsigned char c; 
    in >> c;
    imageImported.tabImage[i] = c;    
  }  
  return imageImported;
}




void 
exportPGM3D(ostream &out, ImagePGM3D image){
    
  out << "P3d" << endl;
    //      << "#Creat by pgm2pgm3d  " 
    //  << "(kerautret@loria.fr)" << endl;
  out << image.width << " " << image.height << " " << image.nbSlice <<endl
      << image.maxValue << endl;  
  
  for(int i=0; i<image.width*image.height*image.nbSlice; i++){
    out << (unsigned char)(image.tabImage[i]);    
  }
  out << endl;
}


