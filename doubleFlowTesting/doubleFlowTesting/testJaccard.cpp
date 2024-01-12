#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>  
#include <string>
#include <climits>
#include <chrono>
#include <opencv2/opencv.hpp>

using namespace std::chrono;
using namespace std;
using namespace cv;


double firstMetric(Mat* frame1, Mat* frame2)
{
    unsigned int correct = 0;

    for (int i = 0; i < frame1->cols; i++)
    {
        for (int j = 0; j < frame1->rows; j++)
        {
            // for a single channel grey scale image
            Scalar intensity1 = frame1->at<uchar>(j, i);
            Scalar intensity2 = frame2->at<uchar>(j, i);
            if (intensity1.val[0] == intensity2.val[0]) correct += 1;
        }
    }
    //cout << "correct: " << correct << endl;
    //cout << "all: " << frame1->cols * frame1->rows << endl;
    return correct/double(frame1->cols * frame1->rows);
}

double jaccardMetric(Mat* frame1, Mat* frame2, string path) // for object 
{
    unsigned int intersect = 0;
    unsigned int unionsect = 0;

    for (int i = 0; i < frame1->cols; i++)
    {
        for (int j = 0; j < frame1->rows; j++)
        {
            // for a single channel grey scale image
            Scalar intensity1 = frame1->at<uchar>(j, i);
            Scalar intensity2 = frame2->at<uchar>(j, i);
            if (intensity1.val[0] == intensity2.val[0] && int(intensity2.val[0])==255) intersect += 1;
            if (int(intensity1.val[0]) == 255 || int(intensity2.val[0]) == 255) unionsect += 1;
        }
    }

    //cout << "\nintersect: " << intersect << endl;
    //cout << "unionsect: " << unionsect << endl;
    //cout << "all:       " << frame1->cols * frame1->rows << endl;
    return intersect / double(unionsect);
}


int main()
{

    vector <string> image_pathes_s1l0 = { "banana3" };
    vector <string> image_pathes_s2l0 = { "banana1", "banana2",  "book", "bool", "bush", 
        "ceramic", "cross", "doll", "elefant", "flower", "fullmoon", "grave", "llama", 
        "memorial", "music", "person1", "person2", "person3", "person4", "person5", 
        "person6", "person7", "person8", "scissors", "sheep", "stone1", "stone2", "teddy", "tennis", 
        "banana3" };

    vector <string> jmZero1 = { "bool",  "fullmoon",  "grave",  "person2",  "person3",  "person4",  "person7", "person8",  "sheep",  "stone1" };
    string bmpSL = "_s2l0.bmp";
    for (string image_path : image_pathes_s2l0)
    {
        if (image_path == "banana3")
            bmpSL = "_s1l0.bmp";
        Mat result = imread("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\results\\" + image_path + bmpSL, IMREAD_GRAYSCALE);
        Mat ideal = imread("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\results\\" + image_path + ".bmp", IMREAD_GRAYSCALE);
        int width = result.size().width;
        int height = result.size().height;
        //cout << m_imageWidth << " " << m_imageHeight << endl;

        // Check for failure
        if (result.empty() || ideal.empty())
        {
            cout << "Could not open or find the image" << endl;
            cin.get(); //wait for any key press
            return -1;
        }

        cout << image_path << ", ";
        cout << setprecision(3) << firstMetric(&result, &ideal) << ", "; // "1 metric = " << 
        cout << setprecision(3) << jaccardMetric(&result, &ideal, image_path) << endl; // << "Jaccard's metric = "
    }

    return 0;
}