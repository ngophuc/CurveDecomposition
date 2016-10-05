#include <iostream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "myfunctions.h"
#include "testfunctions.h"

namespace po=boost::program_options;
using namespace std;
const std::string version="1.0.1";
using namespace std;
using namespace DGtal;

int main(int argc, char *argv[])
{
    /* Read program parameters */
    double angleTol=10.0;
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
            ("help,h", "display this message")
            ("input,i", po::value<std::string>(), "input filename")
            ("imageneDir,d", po::value<std::string>(), "specify the ImaGene directory")
            ("sourceImageWidth", po::value<unsigned int>(), "source image width")
            ("sourceImageHeight", po::value<unsigned int>(), "source image height")
            ("output,o", po::value<std::string>()->default_value("./"), "output filename")
            ("version,v", "display the version number")
            ("eps,e","set output with eps format (default svg)")
            ("maxScale",  po::value<double>()->default_value(15.0), "set the maximal thickness for meaningful thickness detection")
            ("samplingStep",  po::value<double>()->default_value(1.0), "set the threshold of admissible angle in the tangent space")
            ("alphaMax,a",  po::value<double>()->default_value(M_PI/4.0), "set the threshold of admissible angle in the tangent space")
            ("thickness,t",  po::value<double>()->default_value(0.2), "set the tolerance for collinear test in the tangent space")
            ("nbPointCircle,n",  po::value<int>()->default_value(3), "set the minimum number of points on a circle")
            ("isseTol,s",  po::value<double>()->default_value(4.0), "set the ISSE tolerance approximated between segments and arc");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }
    if(vm.count("version")){
        trace.info() << "version: " << version << std::endl;
        return 0;
    }
    po::notify(vm);
    if(vm.count("help")||argc<=1|| !parseOK || !vm.count("input") || !vm.count("output") || !vm.count("imageneDir"))
    {
        trace.info()<< "Contour decomposition." <<std::endl << "Options: "<<std::endl
                    << general_opt << "\n";
        return 0;
    }
    string ImaGeneDIR = vm["imageneDir"].as<std::string>();
    double maxMT = vm["maxScale"].as<double>();
    double stepMT = vm["samplingStep"].as<double>();
    string inputFile=vm["input"].as<std::string>();
    string outputFile=vm["output"].as<std::string>();
    double alphaMax=vm["alphaMax"].as<double>();
    double thicknessMP=vm["thickness"].as<double>();
    double iseTol=vm["isseTol"].as<double>();
    int nbPointCircle=vm["nbPointCircle"].as<int>()-1;
    bool displayImageCanvas=vm.count("sourceImageWidth") && vm.count("sourceImageHeight");
    bool eps = vm.count("eps");
    unsigned int widthCanvas, heightCanvas=0;
    if (displayImageCanvas){
        widthCanvas=vm["sourceImageWidth"].as<unsigned int>();
        heightCanvas=vm["sourceImageHeight"].as<unsigned int>();
    }
    /* Read program parameters */

    /* Read contour */
    std::vector<std::vector<DGtal::Z2i::RealPoint> > vecTmp=PointListReader<DGtal::Z2i::RealPoint>::getPolygonsFromFile(inputFile);
    std::vector<std::vector<DGtal::Z2i::RealPoint> > vecPts=getContours(vecTmp);

    string outVecPts=outputFile+(eps? "OutPts.eps": "OutPts.svg");
    Board2D aBoard;
    aBoard.setLineWidth(100.0);
    aBoard.setFillColor(DGtal::Color::None);
    aBoard.setPenColor(DGtal::Color::Blue);
    for(auto it:vecPts)
        for(auto it_bis:it)
            aBoard.drawCircle(it_bis[0],it_bis[1],1);
    if(eps)
        aBoard.saveEPS(outVecPts.c_str());
    else
        aBoard.saveSVG(outVecPts.c_str());
    /* Read contour */

    /* Compute ATC decomposition */
    //call MT prog
    std::stringstream instruction;
    vector<vector<double> > vecMT;
    for (size_t it_contour = 0;it_contour != vecPts.size();it_contour++)
    {
        std::string contourFile=outputFile+"_"+std::to_string(it_contour)+".sdp";
        writeFile(vecPts.at(it_contour),contourFile.c_str());
        std::string noiseLevelMTFile=outputFile+"MeanThickness_"+std::to_string(it_contour)+".txt";
        instruction << ImaGeneDIR << "/build/tests/TestCompNoiseDetect/displayNoiseBS -srcPolygon " << contourFile
                    << " 0 1 CLOSED -setSampling "<<maxMT<<" "<<stepMT
                    << " -exportNoiseLevel "<< noiseLevelMTFile.c_str();
        std::system(instruction.str().c_str());
        vector<double> oneVecMT=readMeanindfulThicknessFile(noiseLevelMTFile.c_str());
        vecMT.push_back(oneVecMT);
    }
    // build ATC
    std::string outATC=outputFile+(eps? "ATC.eps": "ATC.svg");
    vector<vector<AlphaThickSegmentComputer2DD> > tangentCorverSet=adaptiveTangentCoverDecomposition(vecPts,vecMT,outATC.c_str());
    /* Compute ATC decomposition */

    /* Detect dominant points + optimization */
    vector<vector<RealPoint> > DP;
    std::string outDP=outputFile+(eps? "DP.eps": "DP.svg");
    DP=dominantPointDetection(tangentCorverSet,vecPts,outDP.c_str(),widthCanvas,heightCanvas);
    //Find index of DP on the curve
    vector<vector<int> > indexDP;
    for(size_t it_contour=0; it_contour<vecPts.size(); it_contour++)
    {
        vector<int> index;
        for(vector<RealPoint>::const_iterator it=DP.at(it_contour).begin(); it != DP.at(it_contour).end(); it++)
            index.push_back(findElement(vecPts.at(it_contour),*it));
        indexDP.push_back(index);
    }
    vector<vector<RealPoint> > newDP;
    std::string outNewDP=outputFile+(eps? "newDP.eps": "newDP.svg");
    newDP=dominantPointSimplification(DP,indexDP,vecPts,outNewDP);
    vector<vector<int> > indexNewDP;
    for(size_t it_contour=0; it_contour<vecPts.size(); it_contour++)
    {
        vector<int> index;
        for(vector<RealPoint>::const_iterator it=newDP.at(it_contour).begin(); it != newDP.at(it_contour).end(); it++)
            index.push_back(findElement(vecPts.at(it_contour),*it));
        indexNewDP.push_back(index);
    }
    /* Detect dominant points + optimization */

    /* transform DP into the tangent space */
    vector<vector<RealPoint> > MP=tangentSpaceTransform(newDP,vecPts);
    /* transform DP into the tangent space */

    /* Process arc and seg decomposition */
    vector<vector<double> > DecSegmentsPh,DecArcsPh;
    std::string outDec=outputFile+(eps? "ArcSeg.eps": "ArcSeg.svg");
    std::string outOnlyDec=outputFile+(eps? "OnlyArcSeg.eps": "OnlyArcSeg.svg");
    vector<vector<int> > isolatedPoint=arcSegmentDecomposition(vecPts,indexNewDP,MP,alphaMax,thicknessMP,iseTol,angleTol,nbPointCircle,DecSegmentsPh,DecArcsPh);
    drawDecomposition(vecPts,newDP,isolatedPoint,DecSegmentsPh,DecArcsPh,outDec.c_str(),outOnlyDec.c_str(),widthCanvas,heightCanvas);
    /* Process arc and seg decomposition */

    return 0;
}
