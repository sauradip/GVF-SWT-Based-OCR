# GVF-SWT-Based-OCR
This is a Self Designed OCR , built implementing Gradient Vector Flow and Stroke Width Transform which can Identify any Marathon / Sports Bib Number 

Here main file is gvf_own.m 

I have implemented the Gradient Vector Angle concept to understand the flow of Textual Characters and Then implemented Stroke Width Transform to grow the character from 1 pixel atleast utilising gradient angle 

For missing character , I am finding Orientation of components by Principal Component Analysis and then then moving left or right until i reach another with same orientation 

My code gives around 80 % accuracy for images without tilted text , blurry text , lossy image , any hand/object over the text

Sample Result : 

Input Image :

   ![alt text](https://github.com/sauradip/GVF-SWT-Based-OCR/blob/master/img/t1.jpg)

Detected by my OCR :

![alt text](https://github.com/sauradip/GVF-SWT-Based-OCR/blob/master/img/frame10.jpg)


This result Works best for HD Quality Image but takes Time , whereas this also works best for Low Quality Images as well ( 200x200 ) .
For low Quality Image it performs better than Google Cloud Vision api
