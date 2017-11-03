# GVF-SWT-Based-OCR
This is a Self Designed OCR , built implementing Gradient Vector Flow and Stroke Width Transform which can Identify any Marathon / Sports Bib Number 

Here main file is gvf_own.m 

I have implemented the Gradient Vector Angle concept to understand the flow of Textual Characters and Then implemented Stroke Width Transform to grow the character from 1 pixel atleast utilising gradient angle 

For missing character , I am finding Orientation of components by Principal Component Analysis and then then moving left or right until i reach another with same orientation 

My code gives around 80 % accuracy for images without tilted text , blurry text , lossy image , any hand/object over the text

