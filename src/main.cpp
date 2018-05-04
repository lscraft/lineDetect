#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include "CImg.h"
using namespace std;
using namespace cimg_library;

struct position {
    int row;   // distance
    int col;   // theta
    int num;
};
struct line {
    double k;
    double b;
    int exist;
};
bool cmp(position& a, position& b) {
    return (a.num > b.num);
}
const double pi = 3.14159265359;
const int line_threshold = 700; // for hough space voting
const int edge_threshold = 80; //for edge detector
const int angle = 180;
int mp;

//sobel operator
int filter_Gx[3][3] = { { -1, 0, 1 },{ -2, 0, 2 },{ -1, 0, 1 } };
int filter_Gy[3][3] = { { -1, -2, -1 },{ 0, 0, 0 },{ 1, 2, 1 } };


void sobel(CImg<unsigned char>& pre, CImg<unsigned char>& post) {
    int height = pre.height();
    int width = pre.width();
    post.assign(width, height, 1, 1, 0);
    for (int i = 1; i < width - 1; i++) {
        for (int j = 1; j < height - 1; j++) {
           int temp_hor = 0;
           int temp_ver = 0;
           for (int m = i - 1, a = 0; m <= i + 1; m++, a++) {
                for (int n = j - 1, b = 0; n <= j + 1; n++, b++) {
                    temp_hor += pre(m, n, 0) * filter_Gx[a][b];
                    temp_ver += pre(m, n, 0) * filter_Gy[a][b];
                }
            }
            int temp = sqrt(temp_hor * temp_hor + temp_ver * temp_ver);
            if (temp > edge_threshold) post(i, j, 0) = 255;
            else post(i, j, 0) = 0;
        }
    }
}

//voting
vector<position> vote(CImg<unsigned char>& edge, int** arrvote) {
    double sin_value[angle];
    double cos_value[angle];
    for (int i = 0; i < angle; i++) {
        double di = (double)i;
        sin_value[i] = sin((double)(di / angle) * pi);
        cos_value[i] = cos((double)(di / angle) * pi);
    }
    int height = edge.height();
    int width = edge.width();
    cimg_forXY(edge, i, j) {
        if (edge(i, j) > 0) {
            for (int m = 0; m < angle; m++) {
                double dr = ((double)i) * cos_value[m] + (double)j * sin_value[m];
                int p = round(dr);
                arrvote[m][p + (mp / 2)]++;
            }
        }
    }

    vector<position> pos;
    for (int i = 0; i < angle; i++) {
        for (int j = 0; j < mp; j++) {
           if (arrvote[i][j] > line_threshold) {
                bool newLine = true;
                for (int k = 0; k < pos.size(); k++) {
                    int ok = pos[k].row - mp / 2;
                    //remove close lines
                    if ((abs(pos[k].col - i) < 6 && (abs(pos[k].row - j) < 200) ||  (abs(abs(ok)-abs(j-mp/2)) < 200 && abs(pos[k].col + i-angle) < 6))) {
                        if (pos[k].num < arrvote[i][j]) {
                            pos[k].row = j;
                            pos[k].col = i;
                            pos[k].num = arrvote[i][j];
                        }
                        newLine = false;
                    }

                }
                if (newLine) {
                    position po{ j,i,arrvote[i][j] };
                    pos.push_back(po);
                }
            }
        }
    }
    return pos;
}


//sort vector
void regre(vector<position>& pos) {
    sort(pos.begin(), pos.end(), cmp);
    while (pos.size() > 4) {
        pos.pop_back();
    }
}



vector<line> drawLine(CImg<unsigned char>& drew, vector<position>& pos) {
    vector<line> lines;
    for (int i = 0; i < pos.size(); i++) {
        int x1, x2, y1, y2;
        double ang = (double)(pos[i].col) / angle * pi;
        double si = sin(ang);
        double co = cos(ang);
        if (pos[i].col >= 45 && pos[i].col <= 135) {  //sin bigger than cos
            x1 = 0;
            y1 = (pos[i].row - (mp / 2))/si;
            x2 = drew.width() - 1;
            y2 = ((pos[i].row - (mp / 2)) - (double)x2*co) / si;
        }
        else {
            y1 = 0;
            x1 = (pos[i].row - (mp / 2)) / co;
            y2 = drew.height() - 1;
            x2 = ((pos[i].row - (mp / 2)) - (double)y2*si) / co;
        }
        double k, b;
        k = ((double)y1 - (double)y2) / ((double)x1 - (double)x2);
        b = (double)y1 - k*((double)x1);
        // if K doesn't exist
        if (abs((double)x1 - (double)x2) < 0.001) {
            if (x1 != 0) {
                lines.push_back(line{ 0,(double)x1,0 });
                cout << "line function is x = " << x1 << endl;
            }
            else {
                lines.push_back(line{ 0,(double)x2,0 });
                cout << "line function is x = " << x2 << endl;
            }
        }
        else {
            lines.push_back(line{ k,b,1 });
            cout << "line function is y = " << k << "*x + " << b << endl;
        }
        unsigned char green[] = { 0,255,0 };
        drew.draw_line(x1, y1, x2, y2, green);
    }
    return lines;
}



void drawPoint(CImg<unsigned char>& drew, vector<line>& lines) {
    for (int i = 0; i < lines.size(); i++) {
        for (int j = i + 1; j < lines.size(); j++) {
            double k1 = lines[i].k, b1 = lines[i].b;
            double k2 = lines[j].k, b2 = lines[j].b;
            int x = (b2 - b1) / (k1 - k2);
            int y = k1*(double)x + b1;
            if (lines[i].exist != lines[j].exist && lines[i].exist == 0) {
                x = b1;
                y = k2*(double)x + b2;
            }
            if (lines[i].exist != lines[j].exist && lines[i].exist != 0) {
                x = b2;
                y = k1*(double)x + b1;
            }
            if (x < drew.width() && x >= 0 && y < drew.height() && y >= 0) {
                unsigned char red[] = { 255,0,0 };
                drew.draw_circle(x, y, 30, red);
            }
        
        }
    }
}
void start(const char* filename) {
    CImg<unsigned char> src(filename);  //filename
    CImg<unsigned char> gray = src.get_RGBtoYCbCr().get_channel(0); // to gray scale

    CImg<unsigned char> post_filtrer;



    post_filtrer = gray;
    post_filtrer.blur(2);  // Gaussian filter


    CImg<unsigned char> edge;

    sobel(post_filtrer, edge);       //sobel, get edge of image


    mp = sqrt(2.0)*(double)(edge.width()> edge.height() ? edge.width() : edge.height());
    mp = 2 * mp;                                // the width of hough space

    int* arrvote[angle];

    for (int i = 0; i < angle; i++) {
        arrvote[i] = new int[mp];
        for (int j = 0; j < mp; j++) arrvote[i][j] = 0;
    }


    vector<position> pos = vote(edge, arrvote);  // voting in hough space
    
    for (int i = 0; i < angle; i++) {            //space release;
        delete[]arrvote[i];
    }

    regre(pos);                                  // remove extra positions

    vector<line> lines = drawLine(src, pos);      // transform back to normal space

    drawPoint(src, lines);                     //find the cross point of lines          

    src.display("after");

}
int main() {
    start("test.bmp");
    return 0;
}
