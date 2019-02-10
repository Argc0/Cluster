# Clustering    
This program is an implementation of **K-means** / **K-medoids** **Clustering** algorithms in language **_C++_**. The program implements algorithms for vector clustering in d-dimensional space, using 12 combinations of the following variations. **_Euclidean_** and **_Cosine distances_** are used.      


**_Initialization_**  
1.Random selection of k points (simplest)  
2.K-means++  
**_Assignment_**  
1.Lloyd’s assignment  
2.Assignment by Range search with LSH  
3.Assignment by Range search with Hypercube  
**_Update_**  
1.K-means  
2.Partitioning Around Medoids (PAM) improved like Lloyd’s    

**_How to run(using terminal):_**  
./cluster -i \<input file> -c \<configuration file> -o \<output file> -d \<metric>
