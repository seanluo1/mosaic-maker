/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"

using namespace std;

Point<3> convertToLAB(HSLAPixel pixel) {
    Point<3> result(pixel.h/360, pixel.s, pixel.l);
    return result;
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */

    // vector that stores HSL color values for every pixel
    vector<Point<3>> colors;
    // map that stores unique tile images for point-keys
    map<Point<3>, TileImage*> pixelToTile;

    // convert all images into single HSL points and place them into colors vector
    // also inserts point-image pair into map
    for (unsigned i = 0; i < theTiles.size(); i++) {
        Point<3> p = convertToLAB(theTiles.at(i).getAverageColor());
        colors.push_back(p);
        pixelToTile.insert(make_pair(p, &theTiles.at(i)));
    }
    
    // make a KDTree to represent our image-points
    KDTree<3> myTree = KDTree<3>(colors);
    // make a blank mosaic with the dimensions of source image
    // do i need to delete this later?
    MosaicCanvas* myMosaic = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

    // fill up mosaic with images that closest match the region color of source image
    for (int i = 0; i < theSource.getRows(); i++) {
        for (int j = 0; j < theSource.getColumns(); j++) {
            HSLAPixel regionColor = theSource.getRegionColor(i, j);
            Point<3> p = convertToLAB(regionColor);
            Point<3> nearestNeighbor = myTree.findNearestNeighbor(p);
            myMosaic->setTile(i, j, pixelToTile.find(nearestNeighbor)->second); // retrieve image from map and draw onto mosaic
        }
    }
    return myMosaic;
}

TileImage* get_match_at_idx(const KDTree<3>& tree,
                                  map<Point<3>, int> tile_avg_map,
                                  vector<TileImage>& theTiles,
                                  const SourceImage& theSource, int row,
                                  int col)
{
    // Create a tile which accurately represents the source region we'll be
    // using
    HSLAPixel avg = theSource.getRegionColor(row, col);
    Point<3> avgPoint = convertToLAB(avg);
    Point<3> nearestPoint = tree.findNearestNeighbor(avgPoint);

    // Check to ensure the point exists in the map
    map< Point<3>, int >::iterator it = tile_avg_map.find(nearestPoint);
    if (it == tile_avg_map.end())
        cerr << "Didn't find " << avgPoint << " / " << nearestPoint << endl;

    // Find the index
    int index = tile_avg_map[nearestPoint];
    return &theTiles[index];

}
