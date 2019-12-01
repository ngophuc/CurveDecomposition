#include "decompositionAlgorithm.h"

void drawCanvas(Board2D &aBoard, unsigned int w, unsigned int h){
    aBoard.setLineWidth(0);
    aBoard.drawLine(0, 0, w, 0);
    aBoard.drawLine(w, 0, w, h);
    aBoard.drawLine(w, h, 0, h);
    aBoard.drawLine(0, h, 0, 0);
}

/*************************************/
/*** Burred segments decomposition ***/
/*************************************/
vector<AlphaThickSegmentComputer2DD> blurredSegmentDecomposition(const vector<RealPoint>& aContour, double thickness, const char* filename)
{
    vector<AlphaThickSegmentComputer2DD> fuzzySegmentSet;
    //run over the points on the contours
    for(vector<RealPoint>::const_iterator it=aContour.begin();it!=aContour.end();it++)
    {
        AlphaThickSegmentComputer2DD aSegment(thickness);
        aSegment.init(it);
        /* travel over the contour points and add to the seg */
        while (aSegment.end()!=aContour.end() && aSegment.extendFront()){;}
        if(it==aContour.begin())
            fuzzySegmentSet.push_back(aSegment);
        else if(findElement(aContour,getEndPoint(aSegment))>findElement(aContour,getEndPoint(fuzzySegmentSet.back())))
            fuzzySegmentSet.push_back(aSegment);
        if(getEndPoint(aSegment)==aContour.back() || it==aContour.end())
            break;
    }

    if(filename!=NULL)
    {
        std::string n (filename);
        std::string outputExt = n.substr(n.find_last_of(".")+1);
        Board2D aBoard;
        /* display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for(vector<RealPoint>::const_iterator it=aContour.begin(); it!=aContour.end(); it++)
            aBoard << *it;
        /* display the boundary */
        /* Display boundingbox */
        for(vector<AlphaThickSegmentComputer2DD>::const_iterator it=fuzzySegmentSet.begin();it!=fuzzySegmentSet.end();it++)
            aBoard << SetMode((*it).className(), "BoundingBox")
                   <<(*it);
        /* Display boundingbox */
        if(outputExt=="svg"){
            aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
            aBoard.saveEPS(filename);
        }
    }
    return fuzzySegmentSet;
}
/*************************************/
/*** Burred segments decomposition ***/
/*************************************/

/*********************************************************/
/********* Adaptive Tangent Cover computation ************/
/*********************************************************/
vector<AlphaThickSegmentComputer2DD> adaptiveTangentCoverDecomposition(const vector<RealPoint>& aContour, const vector<double>& vecMT)
{
    //1. Find vector of thickness element
    vector<double> meaningThicknessElement;
    meaningThicknessElement.push_back(vecMT.front());
    for(vector<double>::const_iterator it = vecMT.begin()+1; it != vecMT.end(); it++)
    {
        double m = (*it);
        if (std::find(meaningThicknessElement.begin(), meaningThicknessElement.end(),m)==meaningThicknessElement.end())
            meaningThicknessElement.push_back(m);
    }
    std::sort(meaningThicknessElement.begin(),meaningThicknessElement.end(),sort_increase);

    //2. Compute different thickness tangent covers (blurred segments)
    vector<vector<AlphaThickSegmentComputer2DD> > meaningThicknessTangentCover(meaningThicknessElement.size());
    int index = 0;
    for(vector<double>::const_iterator it = meaningThicknessElement.begin(); it != meaningThicknessElement.end(); it++)
    {
        double thickness = (*it)*sqrt(2);
        cout<<"thickness="<<thickness<<endl;
        vector<AlphaThickSegmentComputer2DD> fuzzySegmentSet = blurredSegmentDecomposition(aContour,thickness,NULL);
        cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
            meaningThicknessTangentCover[index].push_back(*it_bis);
        index++;
    }

    //3. Update thickness of points with tangent covers
    vector<double> vecMTmodified;
    for(vector<double>::const_iterator it = vecMT.begin(); it != vecMT.end(); it++)
        vecMTmodified.push_back(*it);
    for(int it=meaningThicknessTangentCover.size()-1; it>=0; it--)
    {
        vector<AlphaThickSegmentComputer2DD> fuzzySegmentSet = meaningThicknessTangentCover.at(it);//*it;
        double thickness = meaningThicknessElement.at(it);
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            int idStart = findElement(aContour,getStartPoint(*it_bis));
            int idEnd = findElement(aContour,getEndPoint(*it_bis),idStart);
            if(idStart != -1 && idEnd != -1)
            {
                double maxThickness = -1;
                for(int i=idStart; i<=idEnd; i++)
                {
                    double thicknessPoint = vecMT.at(i);
                    if(thicknessPoint > maxThickness)
                        maxThickness = thicknessPoint;
                }
                for(int i=idStart; i<=idEnd; i++)
                {
                    if(maxThickness==thickness)//vecMTmodified.at(i) < maxThickness &&
                        vecMTmodified.at(i) = maxThickness;
                }
            }
            else
                cout<<"negatif"<<endl;
        }
    }

    //4. Travel over the tangent covers and select the segments w.r.t the associated thickness of points
    vector<vector<AlphaThickSegmentComputer2DD> > adaptiveMeaningThicknessTangentCover;
    int idCover = 0;
    for(vector<vector<AlphaThickSegmentComputer2DD> >::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        adaptiveMeaningThicknessTangentCover.push_back(vector<AlphaThickSegmentComputer2DD>());
        vector<AlphaThickSegmentComputer2DD> fuzzySegmentSet = *it;
        vector<AlphaThickSegmentComputer2DD> AdaptiveFuzzySegmentSet;
        int idSeg=0;
        double thickness = meaningThicknessElement.at(idCover);
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            int idStart = findElement(aContour,getStartPoint(*it_bis));
            int idEnd = findElement(aContour,getEndPoint(*it_bis));
            if(idStart != -1 && idEnd != -1)
            {
                bool isGoodMTmodif = false, isGoodMT = false;//true
                for(int i=idStart; i<=idEnd; i++)
                {
                    double thicknessMT = vecMT.at(i); //all elt have same meaningful thickness value (dont contain other meaningful thickness)
                    if(thicknessMT == thickness)
                        isGoodMT = true;
                    double thicknessMTmodif = vecMTmodified.at(i);
                    if(thicknessMTmodif == thickness) //there exist at least one elt in modif having meaningful thickness value
                        isGoodMTmodif = true;
                }
                if(isGoodMTmodif && isGoodMT)
                    AdaptiveFuzzySegmentSet.push_back(*it_bis);
            }
            idSeg++;
        }
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it_bis = AdaptiveFuzzySegmentSet.begin();it_bis != AdaptiveFuzzySegmentSet.end();it_bis++)
            adaptiveMeaningThicknessTangentCover[idCover].push_back(*it_bis);
        idCover++;
    }
    for(vector<vector<AlphaThickSegmentComputer2DD> >::reverse_iterator it1 = adaptiveMeaningThicknessTangentCover.rbegin(); it1 != adaptiveMeaningThicknessTangentCover.rend(); ++it1)
    {
        vector<AlphaThickSegmentComputer2DD>& segmentSet1 = *it1;
        for(vector<vector<AlphaThickSegmentComputer2DD> >::reverse_iterator it2 = it1+1; it2 != adaptiveMeaningThicknessTangentCover.rend(); ++it2)
        {
            vector<AlphaThickSegmentComputer2DD>& segmentSet2 = *it2;
            for (vector<AlphaThickSegmentComputer2DD>::iterator itSeg1 = segmentSet1.begin();itSeg1 != segmentSet1.end();itSeg1++)
            {
                int idCurrentStart = findElement(aContour,getStartPoint(*itSeg1));
                int idCurrentEnd = findElement(aContour,getEndPoint(*itSeg1),idCurrentStart);
                for (vector<AlphaThickSegmentComputer2DD>::iterator itSeg2 = segmentSet2.begin();itSeg2 != segmentSet2.end();itSeg2++)
                {
                    int idStart = findElement(aContour,getStartPoint(*itSeg2));
                    int idEnd = findElement(aContour,getEndPoint(*itSeg2),idStart);
                    if(idCurrentStart<=idStart && idCurrentEnd>=idEnd)
                    {
                        segmentSet2.erase(itSeg2);
                        itSeg2--;
                    }
                }
            }
        }
    }

    //5. Reorder the multi-thickness tangent cover
    vector<AlphaThickSegmentComputer2DD> adaptiveTangentCover;
    int seg=0,nbSeg=0;
    vector<int> idThicknessCover; //stock idSeg of the last seg at idThicknessCover
    for(int it=0; it<meaningThicknessElement.size();it++)
        idThicknessCover.push_back(0);
    for(int it = 0; it < adaptiveMeaningThicknessTangentCover.size(); it++)
        nbSeg += (adaptiveMeaningThicknessTangentCover.at(it)).size();

    while (seg<nbSeg)
    {
        int idMinStart = aContour.size(), idMinEnd = aContour.size();
        int idMin=-1, idSeg=-1;
        //scan adaptiveMeaningThicknessTangentCover
        for(int it = 0; it < adaptiveMeaningThicknessTangentCover.size(); it++)//thickness level = it
        {
            //current seg of thickness level idThicknessCover.at(i)
            int idCurrentSeg = idThicknessCover.at(it);
            if(idCurrentSeg<(adaptiveMeaningThicknessTangentCover.at(it)).size())
            {
                //get idStart and idEnd of seg
                int idCurrentStart = findElement(aContour,getStartPoint((adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg)));
                int idCurrentEnd = findElement(aContour,getEndPoint((adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg)),idCurrentStart);
                if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size())
                {
                    //find min idCurrentStart
                    if(idMinStart==idCurrentStart && idMinEnd<idCurrentEnd)
                    {
                        if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size()-1)
                        {
                            idThicknessCover.at(idMin) = idThicknessCover.at(idMin) + 1;
                            seg++;
                        }
                        idSeg = idCurrentSeg;
                        idMin = it;
                        idMinStart = idCurrentStart;
                        idMinEnd = idCurrentEnd;
                    }
                    else if(idMinStart>idCurrentStart && idMinEnd>=idCurrentEnd)
                    {
                        idSeg = idCurrentSeg;
                        idMin = it;
                        idMinStart = idCurrentStart;
                        idMinEnd = idCurrentEnd;
                    }
                }
            }
        }
        adaptiveTangentCover.push_back((adaptiveMeaningThicknessTangentCover.at(idMin)).at(idSeg));
        idThicknessCover.at(idMin) = idThicknessCover.at(idMin) + 1;
        seg++;
    }

    return adaptiveTangentCover;
}

vector<vector<AlphaThickSegmentComputer2DD> > adaptiveTangentCoverDecomposition(const vector<vector<RealPoint> >& aContour, const vector<vector<double> >& vecMT, std::string filename, unsigned int w, unsigned int h)
{
    vector<vector<AlphaThickSegmentComputer2DD> > vecATC;
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        vecATC.push_back(adaptiveTangentCoverDecomposition(aContour.at(it_contour),vecMT.at(it_contour)));

    Board2D aBoard;
    if(w!=0 && h!=0)
        drawCanvas(aBoard, w, h);
    GradientColorMap<int> cmap_grad( 0, aContour.size()+1 );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );
    cmap_grad.addColor( Color( 25, 255, 255 ) );
    cmap_grad.addColor( Color( 255, 25, 255 ) );
    cmap_grad.addColor( Color( 25, 25, 25 ) );
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it=vecATC.at(it_contour).begin(); it != vecATC.at(it_contour).end(); ++it)
            aBoard << SetMode((*it).className(), "BoundingBox")
                   << CustomStyle("AlphaThickSegment/BoundingBox",  new CustomPenColor(cmap_grad(it_contour)))
                   << *it;
    std::string outputExt = filename.substr(filename.find_last_of(".")+1);
    if(outputExt=="svg")
        aBoard.saveSVG(filename.c_str());
    else if (outputExt == "eps")
        aBoard.saveEPS(filename.c_str());

    return vecATC;
}
/*********************************************************/
/********* Adaptive Tangent Cover computation ************/
/*********************************************************/

/************************************/
/**** Dominant points detections ****/
/************************************/
vector<RealPoint> dominantPointDetection(const vector<AlphaThickSegmentComputer2DD>& vecTangentCover, const vector<RealPoint>& aContour)
{
    int p=1, q=0, Eq=0, Bp=0, m=vecTangentCover.size();
    vector<RealPoint> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        Bp = findStartElement(aContour,vecTangentCover[p]);
        Eq = findEndElement(aContour,vecTangentCover[q]);
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findStartElement(aContour,vecTangentCover[p]);
            Eq = findEndElement(aContour,vecTangentCover[q]);
            p++;
        }
        if(Eq<Bp)
        {
            p--;//when out of loop => Eq<Bp
            pile.push_back(RealPoint(q,p-1));
            q=p-1;
        }
        else
        {
            if(p>=m || q>=m)
                pile.push_back(RealPoint(q,p-1));
            else
                q = p;
        }
    }
    /* pile construction */

    /* calcul point dominant */
    vector<RealPoint> DP;
    for(vector<RealPoint>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        //find point with max (local) curvature
        int indexB = findStartElement(aContour,vecTangentCover[p]);
        int indexE = findEndElement(aContour,vecTangentCover[q]);
        int indexC = indexB;
        int indexBi = findStartElement(aContour,vecTangentCover[q]);
        int indexEi = findEndElement(aContour,vecTangentCover[p]);

        for(int i = indexB+1; i<=indexE; i++)
        {
            int modIndex = i;
            if(fabs(acuteAngle(aContour.at(indexBi),aContour.at(modIndex),aContour.at(indexEi))) <=
                    fabs(acuteAngle(aContour.at(indexBi),aContour.at(indexC),aContour.at(indexEi))))
                indexC = modIndex;
        }

        if(DP.empty() || aContour[indexC] != DP.back())
            DP.push_back(aContour[indexC]);
    }
    DP.insert(DP.begin(),aContour.front());
    DP.insert(DP.end(),aContour.back());
    /* calcul point dominant */

    return DP;
}

vector<vector<RealPoint> > dominantPointDetection(const vector<vector<AlphaThickSegmentComputer2DD> >& vecTangentCover, const vector<vector<RealPoint> >& aContour, string filename, unsigned int w, unsigned int h)
{
    vector<vector<RealPoint> > vecDP;
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        //cout<<"contour "<<it_contour<<": size="<<aContour.at(it_contour).size()<<endl;
        vecDP.push_back(dominantPointDetection(vecTangentCover.at(it_contour),aContour.at(it_contour)));

    Board2D aBoard;
    if(w!=0 && h!=0)
        drawCanvas(aBoard, w, h);
    aBoard << SetMode("PointVector", "Grid");
    GradientColorMap<int> cmap_grad( 0, aContour.size()+1 );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );
    cmap_grad.addColor( Color( 25, 255, 255 ) );
    cmap_grad.addColor( Color( 255, 25, 255 ) );
    cmap_grad.addColor( Color( 25, 25, 25 ) );
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
    {
        /* Display the boundary */
        for (vector<RealPoint>::const_iterator it = aContour.at(it_contour).begin(); it != aContour.at(it_contour).end(); it++)
            aBoard << *it;
        /* Display the boundary */
        /* Display the segments by DP */
        aBoard.setPenColor(cmap_grad(it_contour));
        aBoard.setLineWidth(100.0);
        for(vector<RealPoint>::const_iterator it = vecDP.at(it_contour).begin(); it+1 != vecDP.at(it_contour).end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        /* Display the segments by DP */
    }

    std::string outputExt = filename.substr(filename.find_last_of(".")+1);
    if(outputExt=="svg")
        aBoard.saveSVG(filename.c_str());
    else if (outputExt == "eps")
        aBoard.saveEPS(filename.c_str());

    return vecDP;
}
/************************************/
/**** Dominant points detections ****/
/************************************/

/***********************************/
/**** Dominant points selection ****/
/***********************************/
vector<RealPoint> dominantPointSimplification(const vector<RealPoint>& DP, const vector<int>& indexDP, const vector<RealPoint>& aContour, int nbMinPts)
{
    if(DP.size()<=nbMinPts)
        return DP;
    vector<RealPoint> selectedDP;
    for(vector<RealPoint>::const_iterator it = DP.begin(); it != DP.end(); it++)
        selectedDP.push_back(*it);
    vector<int> indexSelectedDP;
    for(vector<int>::const_iterator it = indexDP.begin(); it != indexDP.end(); it++)
        indexSelectedDP.push_back(*it);
    vector<double> iseDP;
    double ise = error_ISE(aContour,selectedDP,indexSelectedDP,false);
    double cr = error_CR(aContour,selectedDP);
    double error_prev, error;
    error_prev = error = (cr*cr)/ise;
    int lastIndex=-1,lastIndexSelectedDP=-1;
    RealPoint lastSelectedDP;

    while(error_prev <= error){
        error_prev = error;
        iseDP.push_back(aContour.size()*aContour.size());//Keep the fist point as DP
        //cout<<"iseDP first="<<iseDP.back()<<endl;
        for(int it = 1; it+1 < selectedDP.size(); it++)
        {
            double ise = 0.0, d = 0.0, angle = 0.0;
            int indexStart = indexSelectedDP.at(it-1);
            int indexEnd = indexSelectedDP.at(it+1);
            if(indexStart < 0 || indexEnd < 0)
                cout<<"Pb in error_ISE index of DP"<<endl;
            for(int index = indexStart+1; index < indexEnd; index++)
            {
                d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.at(indexEnd));
                ise += d*d;
            }
            angle = acuteAngle(aContour.at(indexStart),
                               aContour.at(indexSelectedDP.at(it)),
                               aContour.at(indexEnd));
            iseDP.push_back(ise/angle);
            //cout<<"it "<<it<<" : iseDP "<<iseDP.back()<<endl;
        }
        iseDP.push_back(aContour.size()*aContour.size());//Keep the fist point as DP
        //cout<<"iseDP last="<<iseDP.back()<<endl;

        vector<int> indexOrderedDP = absSortIndex(iseDP,true);
        /* erease the first elt */
        lastIndex = indexOrderedDP.front();
        lastSelectedDP = selectedDP.at(lastIndex);
        lastIndexSelectedDP = indexSelectedDP.at(lastIndex);
        //cout<<"=> min ise at "<<indexOrderedDP.front()<<endl;
        selectedDP.erase(selectedDP.begin()+indexOrderedDP.front());
        indexSelectedDP.erase(indexSelectedDP.begin()+indexOrderedDP.front());
        iseDP.clear();
        /* erease the first elt */

        /* update errors */
        ise = error_ISE(aContour,selectedDP,indexSelectedDP,false);
        cr = error_CR(aContour,selectedDP);
        error = (cr*cr)/ise;
        /* update errors */
    }
    /* Push back the last element */
    if(lastIndex !=-1)
    {
        selectedDP.insert(selectedDP.begin()+lastIndex,lastSelectedDP);
        indexSelectedDP.insert(indexSelectedDP.begin()+lastIndex,lastIndexSelectedDP);
    }
    /* Push back the last element */

    return selectedDP;
}

vector<vector<RealPoint> > dominantPointSimplification(const vector<vector<RealPoint> > &DP, const vector<vector<int> >& indexDP, const vector<vector<RealPoint> > &aContour, string filename, int nbMinPts, unsigned int w, unsigned int h)
{
    vector<vector<RealPoint> > selectedDP;
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        selectedDP.push_back(dominantPointSimplification(DP.at(it_contour),indexDP.at(it_contour),aContour.at(it_contour),nbMinPts));

    Board2D aBoard;
    if(w!=0 && h!=0)
        drawCanvas(aBoard, w, h);
    /* Display the boundary */
    aBoard << SetMode("PointVector", "Both");
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
    {
        for (vector<RealPoint>::const_iterator it = aContour.at(it_contour).begin(); it != aContour.at(it_contour).end(); it++)
            aBoard << *it;
        /* Display the boundary */
        aBoard.setLineWidth(100.0);
        /* Display the segments by old DP *
        aBoard.setPenColor(Color(192, 0, 0));
        for(vector<RealPoint>::const_iterator it = DP.at(it_contour).begin(); it+1 != DP.at(it_contour).end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        /* Display the segments by old DP */
        /* Display the segments by new DP */
        aBoard.setPenColor(Color(0, 192, 0));
        for(vector<RealPoint>::const_iterator it = selectedDP.at(it_contour).begin(); it+1 != selectedDP.at(it_contour).end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        /* Display the segments by new DP */
    }
    std::string outputExt = filename.substr(filename.find_last_of(".")+1);
    if(outputExt=="svg")
        aBoard.saveSVG(filename.c_str());
    else if (outputExt == "eps")
        aBoard.saveEPS(filename.c_str());

    return selectedDP;
}
/************************************/
/**** Dominant points selection *****/
/************************************/

/**************************************/
/**** Tangent space transformation ****/
/**************************************/
vector<RealPoint> tangentSpaceTransform(const vector<RealPoint>& DP)
{
    vector<RealPoint> MidPoints;
    double l=0.0, a=0.0, totalAngle=0.0;
    vector<double> alpha;
    vector<double> length;
    double totalLength=0;

    //scan the Dominant Points
    alpha.push_back(a);
    for(vector<RealPoint>::const_iterator it=DP.begin(); it+1!=DP.end(); it++)
    {
        //2 consecutive DP=1 segment => store the length and angle between segments
        l=distancePoints(*it,*(it+1));
        length.push_back(l);
        totalLength +=l;
        if(it!=DP.begin())
        {
            RealPoint v1=RealPoint((*it)[0]-(*(it-1))[0], (*it)[1]-(*(it-1))[1]);
            RealPoint v2=RealPoint((*(it+1))[0]-(*it)[0],(*(it+1))[1]-(*it)[1]);
            a=signedAngle(v1,v2);
            alpha.push_back(a);
            totalAngle += a;
        }
    }

    //transform into the tangent space (x,y)=(length,alpha)
    vector<RealPoint> Ti1,Ti2,Ti;
    RealPoint p(0.0,0.0);
    Ti.push_back(p);
    Ti2.push_back(p);

    size_t it;
    for(it=0; it<alpha.size()-1; it++)
    {
        p[0]=Ti2.back()[0]+length[it];
        p[1]=Ti2.back()[1];
        Ti1.push_back(p);
        Ti.push_back(p);
        p[0]=Ti1.back()[0];
        p[1]=Ti1.back()[1]+alpha[it+1];
        Ti2.push_back(p);
        Ti.push_back(p);
    }
    p[0]=Ti2.back()[0]+length[it];
    p[1]=Ti2.back()[1];
    Ti.push_back(p);

    for(vector<RealPoint>::const_iterator it=Ti.begin(); it!=Ti.end(); it+=2)
        MidPoints.push_back(RealPoint(((*it)[0]+(*(it+1))[0])/2.0,((*it)[1]+(*(it+1))[1])/2.0));

    return MidPoints;
}

vector<vector<RealPoint> > tangentSpaceTransform(const vector<vector<RealPoint> > &DP, const vector<vector<RealPoint> > &aContour)
{
    vector<vector<RealPoint> > vecMidPoints;
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        vecMidPoints.push_back(tangentSpaceTransform(DP.at(it_contour)));
    return vecMidPoints;
}
/**************************************/
/**** Tangent space transformation ****/
/**************************************/

/*****************************************************/
/****Decomposition of Curve into Segments and Arcs ***/
/*****************************************************/
int findBestFittingCircle(const vector<RealPoint> aContour, int idStart, int idEnd)
{
    int oneThird=(idEnd-idStart)/3;
    double ise,radius,minIse=-1;
    int idMin=-1;
    RealPoint center;
    for(int idMid=idStart+oneThird ; idMid<(idEnd-oneThird); idMid++)
    {
        center=determineCenter(aContour.at(idStart),aContour.at(idMid%aContour.size()),aContour.at(idEnd%aContour.size()));
        radius=(determineRadius(center,aContour.at(idStart)) + determineRadius(center,aContour.at(idMid%aContour.size())) + determineRadius(center,aContour.at(idEnd%aContour.size())))/3.0;
        ise=iseContourCircle(aContour,idStart,idEnd%aContour.size(),center,radius);
        if((minIse<0) || (ise<minIse))
        {
            minIse=ise;
            idMin=idMid%aContour.size();
        }
    }
    return idMin;
}

//Dealing with iseTol=0/st/20 => only segments/mix segments and arcs/only arcs
vector<int> arcSegmentDecomposition(const vector<RealPoint> &aContour, const vector<int> &indexDP, const vector<RealPoint> &MP, double alphaMax, double thickness, double iseTol, double angleTol, int nbPointCir, vector<double> &segments, vector<double> &arcs)
{
    vector<bool> isArc;
    for(size_t it=0; it<MP.size(); it++)
        isArc.push_back(false);

    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    vector<RealPoint> arcIndex;
    if(MP.size()<2)
        iseTol=0;
    if(iseTol==0.0)//min => only segments => nothing to do
        for(size_t it=0; it<MP.size(); it++)
            isolatedVector.push_back(1);
    else if(iseTol==20.0)//max => only arcs
    {
        for(size_t it=0; it<MP.size(); it++)
            isolatedVector.push_back(0);

        /*********** Decomposition into blurred segments ****/
        vector<RealPoint>::const_iterator it_MP=MP.begin();
        for(size_t it_start=0; it_start<MP.size(); it_start++)
        {
            int it_end=it_start;
            AlphaThickSegmentComputer2DD aSegment(thickness);
            aSegment.init(it_MP);
            while(aSegment.end()!=MP.end() && aSegment.extendFront())
                it_end++;
            it_MP++;
            if(blurredSegmentTS.size() == 0 || (blurredSegmentTS.size() != 0 && (findElement(MP,getEndPoint(aSegment)) > findElement(MP,getEndPoint(blurredSegmentTS.back())))))
            {
                int idEnd,idMid,idBegin=indexDP.at(it_start);
                idEnd=indexDP.at(it_end+1);
                int idEndOld=idEnd;
                if(fabs(MP.at(findElement(MP,getStartPoint(aSegment)))[1]-MP.at(findElement(MP,getEndPoint(aSegment)))[1])>(1.4*M_PI))
                    idEnd=(int)((idBegin+3*idEnd)/4);
                idMid=findBestFittingCircle(aContour,idBegin,idEnd);
                if(idMid==-1)
                    cout<<"idMid==-1 ==> idBegin=="<<idBegin<<" and idEnd="<<idEnd<<" aSegment.getNumberSegmentPoints() "<<aSegment.getNumberSegmentPoints()<<endl;
                RealPoint center=determineCenter(aContour.at(idBegin),
                                                 aContour.at(idMid),
                                                 aContour.at(idEnd));
                double radius=(determineRadius(center,
                                               aContour.at(idBegin)) +
                               determineRadius(center,aContour.at(idMid)) +
                               determineRadius(center,aContour.at(idEnd)))/3.0;
                double ise_Seg=0;
                for(int i=it_start; i<=it_end; i++)
                    ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                arcIndex.push_back(Point(idBegin,idEndOld));
                for(int i=it_start; i<=it_end; i++)
                    isArc[i]=true;

                arcs.push_back(center[0]);
                arcs.push_back(center[1]);
                arcs.push_back(radius);
                double startAngle=atan2(aContour.at(idBegin)[1]-center[1],
                        aContour.at(idBegin)[0]-center[0]);
                double endAngle=atan2(aContour.at(idEndOld)[1]-center[1],
                        aContour.at(idEndOld)[0]-center[0]);
                double neg=isLeft(aContour.at(idBegin),
                                  aContour.at(idMid),
                                  aContour.at(idEnd))<0 ? 1 : -1;
                arcs.push_back(startAngle);
                arcs.push_back(endAngle);
                arcs.push_back(neg);

                blurredSegmentTS.push_back(aSegment);
            }
            if(getEndPoint(aSegment)==MP.back())
                break;
        }
        /*********** Decomposition into blurred segment ****/
    }
    else //in between => mix between arcs and segments
    {
        /*********** Test of isolated points ***********/
        for(size_t it=0; it<MP.size(); it++)
        {
            if(it==0)
            {
                if(fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) //MP[0] is an isolated point => segment !
                    isolatedVector.push_back(1);//SEG
                else//==> arc
                    isolatedVector.push_back(0);//ARC
            }
            else if(it==(int)MP.size()-1)
            {
                if(fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax)//MP[0] is an isolated point => segment !
                    isolatedVector.push_back(1);//SEG
                else // ARC or JUNCTION POINT
                    isolatedVector.push_back(0);//ARC
            }
            else
            {
                if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) &&
                        (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
                    isolatedVector.push_back(1);//SEG
                else
                {
                    if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                        isolatedVector.push_back(-1);//JONCTION
                    else
                        isolatedVector.push_back(0);//ARC
                }
            }
        }
        /*********** Test of isolated points ***********/

        /*********** Decomposition by blurred segment ****/
        vector<RealPoint>::const_iterator it_MP=MP.begin();
        for(size_t it_start=0; it_start<isolatedVector.size(); it_start++)
        {
            int it_end=it_start;
            AlphaThickSegmentComputer2DD aSegment(thickness);
            aSegment.init(it_MP);
            while(it_end<isolatedVector.size() && isolatedVector.at(it_end)!=1 &&
                  it_end<isolatedVector.size()-1 && isolatedVector.at(it_end+1)!=1 && //For the isolated point comes after
                  aSegment.end()!=MP.end() && aSegment.extendFront()) it_end++;
            it_MP++;
            if(aSegment.getNumberSegmentPoints()>=nbPointCir)//at least (nbPointCir+1) points on the circle
            {
                //cout<<"it can be a circle => aSegment.getNumberSegmentPoints()="<<aSegment.getNumberSegmentPoints()<<endl;
                if(blurredSegmentTS.size() == 0 || (blurredSegmentTS.size() != 0 &&
                                                    (findElement(MP,getEndPoint(aSegment)) > findElement(MP,getEndPoint(blurredSegmentTS.back())))))
                {
                    int idEnd,idMid,idBegin=indexDP.at(it_start);
                    idEnd=indexDP.at(it_end+1);//FIXME : it_end
                    int idEndOld=idEnd;
                    if(fabs(MP.at(findElement(MP,getStartPoint(aSegment)))[1]-MP.at(findElement(MP,getEndPoint(aSegment)))[1])>(1.4*M_PI))
                        idEnd=(int)((idBegin+3*idEnd)/4);
                    idMid=findBestFittingCircle(aContour,idBegin,idEnd);
                    if(idMid==-1)//FIXME: idBegin=idEnd in case of closed circle curve !!!
                        cout<<"idMid==-1 ==> idBegin=="<<idBegin<<" and idEnd="<<idEnd<<" aSegment.getNumberSegmentPoints() "<<aSegment.getNumberSegmentPoints()<<endl;
                    double linAngle=relativeAngle(aContour.at(idBegin),
                                                  aContour.at(idMid),
                                                  aContour.at(idEnd))*180/M_PI;
                    RealPoint center=determineCenter(aContour.at(idBegin),
                                                     aContour.at(idMid),
                                                     aContour.at(idEnd));
                    double radius=(determineRadius(center,aContour.at(idBegin)) +
                                   determineRadius(center,aContour.at(idMid)) +
                                   determineRadius(center,aContour.at(idEnd)))/3.0;
                    double ise_Arc=iseContourCircle(aContour,idBegin,idEnd,center,radius);
                    double ise_Seg=0;
                    for(int i=it_start; i<=it_end; i++)//FIXME : i<it_end
                        ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                    if (!(ise_Arc>iseTol*ise_Seg || (180-linAngle)<angleTol))//fabs(180-linAngle)<angleTol
                    {
                        arcIndex.push_back(Point(idBegin,idEndOld));
                        for(int i=it_start; i<=it_end; i++)//FIXME : i<it_end
                            isArc[i]=true;
                        //cout<<"it is a circle"<<endl;
                        arcs.push_back(center[0]);
                        arcs.push_back(center[1]);
                        arcs.push_back(radius);
                        double startAngle=atan2(aContour.at(idBegin)[1]-center[1],
                                aContour.at(idBegin)[0]-center[0]);
                        double endAngle=atan2(aContour.at(idEndOld)[1]-center[1],
                                aContour.at(idEndOld)[0]-center[0]);
                        double neg=isLeft(aContour.at(idBegin),
                                          aContour.at(idMid),
                                          aContour.at(idEnd))<0 ? 1 : -1;
                        arcs.push_back(startAngle);
                        arcs.push_back(endAngle);
                        arcs.push_back(neg);

                        blurredSegmentTS.push_back(aSegment);
                    }
                }
                if(getEndPoint(aSegment)==MP.back())
                    break;
            }
        }
        /*********** Decomposition by blurred segment ****/
    }

    /* Dealing with inclusion */
    if(arcIndex.size()>0)
    {
        vector<int> arcerase;
        for(size_t i=1; i<arcIndex.size()-1; i++)
        {
            int idEndPrev=arcIndex.at(i-1)[1];
            int idBeginSucc=arcIndex.at(i+1)[0];
            if(idEndPrev>=idBeginSucc)
            {
                ///cout<<"remove : "<<i<<" => "<<arcIndex.at(i)[0]<<" , "<<arcIndex.at(i)[1]<<endl;
                arcerase.push_back(i);
                i++;
            }
        }
        for(int i=(int)arcerase.size()-1; i>=0; i--)
        {
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcs.erase(arcs.begin()+(6*arcerase.at(i)));
            arcIndex.erase(arcIndex.begin()+arcerase.at(i));
        }
    }
    /* Dealing with inclusion */

    /* Dealing with intersection */
    if(arcIndex.size()>0)
    {
        vector<int> arcIntersection;
        for(size_t i=0; i<arcIndex.size()-1; i++)
        {
            int idBeginPrev=arcIndex.at(i)[0];
            int idEndPrev=arcIndex.at(i)[1];
            int idBegin=arcIndex.at(i+1)[0];
            int idEnd=arcIndex.at(i+1)[1];
            if(idEndPrev>idBegin)
            {
                arcIntersection.push_back(i);
                int idEndMid=((idEndPrev+idBegin)/2)%aContour.size();
                int idMid1=findBestFittingCircle(aContour,idBeginPrev,idEndMid)%aContour.size();
                RealPoint center1=determineCenter(aContour.at(idBeginPrev),
                                              aContour.at(idMid1),
                                              aContour.at(idEndMid));
                double radius1=(determineRadius(center1,aContour.at(idBeginPrev)) +
                                determineRadius(center1,aContour.at(idMid1)) +
                                determineRadius(center1,aContour.at(idEndMid)))/3.0;
                double startAngle1=atan2(aContour.at(idBeginPrev)[1]-center1[1],
                        aContour.at(idBeginPrev)[0]-center1[0]);
                double endAngle1=atan2(aContour.at(idEndMid)[1]-center1[1],
                        aContour.at(idEndMid)[0]-center1[0]);
                double neg1=isLeft(aContour.at(idBeginPrev),
                                   aContour.at(idMid1),
                                   aContour.at(idEndMid))<0 ? 1 : -1;
                arcs.at(6*i)=center1[0];
                arcs.at(6*i+1)=center1[1];
                arcs.at(6*i+2)=radius1;
                arcs.at(6*i+3)=startAngle1;
                arcs.at(6*i+4)=endAngle1;
                arcs.at(6*i+5)=neg1;
                arcIndex.at(i)[0]=idBeginPrev;
                arcIndex.at(i)[1]=idEndMid;

                int idBeginMid=((idEndPrev+idBegin)/2)%aContour.size();
                int idMid2=findBestFittingCircle(aContour,idBeginMid,idEnd)%aContour.size();
                RealPoint center2=determineCenter(aContour.at(idBeginMid),
                                              aContour.at(idMid2),
                                              aContour.at(idEnd%aContour.size()));
                double radius2=(determineRadius(center2,aContour.at(idBeginMid)) +
                                determineRadius(center2,aContour.at(idMid2)) +
                                determineRadius(center2,aContour.at(idEnd%aContour.size())))/3.0;
                double startAngle2=atan2(aContour.at(idBeginMid)[1]-center2[1],
                        aContour.at(idBeginMid)[0]-center2[0]);
                double endAngle2=atan2(aContour.at(idEnd%aContour.size())[1]-center2[1],
                        aContour.at(idEnd%aContour.size())[0]-center2[0]);
                double neg2=isLeft(aContour.at(idBeginMid),
                                   aContour.at(idMid2),
                                   aContour.at(idEnd%aContour.size()))<0 ? 1 : -1;
                arcs.at(6*(i+1))=center2[0];
                arcs.at(6*(i+1)+1)=center2[1];
                arcs.at(6*(i+1)+2)=radius2;
                arcs.at(6*(i+1)+3)=startAngle2;
                arcs.at(6*(i+1)+4)=endAngle2;
                arcs.at(6*(i+1)+5)=neg2;
                arcIndex.at(i+1)[0]=idBeginMid;
                arcIndex.at(i+1)[1]=idEnd;
            }
        }
    }
    /* Dealing with intersection */

    for(size_t i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
        {
            segments.push_back((aContour.at(indexDP.at(i)))[0]);
            segments.push_back((aContour.at(indexDP.at(i)))[1]);
            segments.push_back((aContour.at(indexDP.at(i+1)))[0]);
            segments.push_back((aContour.at(indexDP.at(i+1)))[1]);
        }
    }

    return isolatedVector;
}

vector<vector<int> > arcSegmentDecomposition(const vector<vector<RealPoint> > &aContour, const vector<vector<int> > &indexDP, const vector<vector<RealPoint> > &MP, double alphaMax, double thickness, double iseTol, double angleTol, int nbPointCir, vector<vector<double> > &segments, vector<vector<double> > &arcs)
{
    vector<vector<int> > vecIsolated;//SEG:1,ARC:0,JOINCTION:-1
    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
    {
        vector<double> DecSegmentsPh,DecArcsPh;
        //cout<<"it_contour "<<it_contour<<": has MP.size()="<<MP.at(it_contour).size()<<endl;
        vector<int> isolated=arcSegmentDecomposition(aContour.at(it_contour),indexDP.at(it_contour),MP.at(it_contour),alphaMax,thickness,iseTol,angleTol,nbPointCir,DecSegmentsPh,DecArcsPh);
        vecIsolated.push_back(isolated);
        segments.push_back(DecSegmentsPh);
        arcs.push_back(DecArcsPh);
    }
    return vecIsolated;
}
/*****************************************************/
/****Decomposition of Curve into Segments and Arcs ***/
/*****************************************************/

/*****************************************************/
/****Draw the decomposition into Segments and Arcs ***/
/*****************************************************/
void drawDecomposition(const vector<RealPoint>& aContour, const vector<RealPoint>& DP, const vector<int> isolatedPoint, const vector<double>& segments, const vector<double>& arcs, Board2D& DecSAboard, Board2D& onlyDecSAboard)
{
    /// Display contour points
    for(vector<RealPoint>::const_iterator it=aContour.begin();it!=aContour.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,0,0), Color(204,204,204)) )
                   << *it;
    /// Display contour points
    /// Display dominant points
    for(vector<RealPoint>::const_iterator it=DP.begin(); it!=DP.end(); it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,0,255), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    //DecSAboard.setPenColor(Color( 0, 0, 0));
    DecSAboard.setPenColor(Color( 0, 255, 0));
    DecSAboard.setLineWidth(50.0);
    onlyDecSAboard.setPenColor(Color( 0, 255, 0));
    onlyDecSAboard.setLineWidth(50.0);
    for(vector<double>::const_iterator it=segments.begin(); it!=segments.end(); it=it+4)
    {
        DecSAboard.drawLine((*it),*(it+1),*(it+2),*(it+3),1);
        onlyDecSAboard.drawLine((*it),*(it+1),*(it+2),*(it+3),1);
    }
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setLineWidth(100.0);
    onlyDecSAboard.setPenColor(Color( 255, 0, 0));
    onlyDecSAboard.setLineWidth(100.0);
    for(vector<double>::const_iterator it=arcs.begin(); it!=arcs.end(); it=it+6)
    {
        Point center=Point((*it),*(it+1));
        double radius=*(it+2);
        double startAngle=*(it+3);
        double endAngle=*(it+4);
        bool neg=(*(it+5))>0 ? true : false;
        DecSAboard.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
        onlyDecSAboard.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
    }
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(size_t it=0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION=ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED=SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
}

void drawDecomposition(const vector<vector<RealPoint> >& aContour, const vector<vector<RealPoint> >& DP, const vector<vector<int> > isolatedPoint, const vector<vector<double> >& segments, const vector<vector<double> >& arcs, string filename, string onlyfilename, unsigned int w, unsigned int h)
{
    Board2D DecSAboard, onlyDecSAboard;
    if(w!=0 && h!=0)
    {
        drawCanvas(DecSAboard, w, h);
        drawCanvas(onlyDecSAboard, w, h);
    }

    for(size_t it_contour=0; it_contour<aContour.size(); it_contour++)
        drawDecomposition(aContour.at(it_contour),DP.at(it_contour),isolatedPoint.at(it_contour),segments.at(it_contour),arcs.at(it_contour),DecSAboard,onlyDecSAboard);

    std::string outputExt = filename.substr(filename.find_last_of(".")+1);
    if(outputExt=="svg")
    {
        DecSAboard.saveSVG(filename.c_str());
        onlyDecSAboard.saveSVG(onlyfilename.c_str());
    }
    else if (outputExt == "eps")
    {
        DecSAboard.saveEPS(filename.c_str());
        onlyDecSAboard.saveEPS(onlyfilename.c_str());
    }
}
/*****************************************************/
/****Draw the decomposition into Segments and Arcs ***/
/*****************************************************/
