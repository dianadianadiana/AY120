import numpy as np

def max_peaks(arr, width, lower_limit):
    """ Returns an array of the indexes where the indexes indicate where the peaks are"""
    #return [i for i in range(1, len(arr)-1) if arr[i-1]<arr[i] and arr[i+1]<arr[i] and arr[i]>lower_limit]
    return [i for i in range(1, len(arr)-1) if all(arr[i] > arr[i-width:i]) and all(arr[i]>arr[i+1:i+width]) and arr[i]>lower_limit] # 36 is arbitrary
    #assuming that there is nothing before 36

def find_peaks(img):

    avg = np.median(img)
    lim = 1.5*avg
    #arr[i][j] == [row][col] ~~~ [y][x] so when we plot them, we need to plot the y values, x values
    # find the initial pixel values that are above the limit (=constant*avg)
    # img_arr has the following syntax:
    # the index of img_arr corresponds to the row value of the pixel
    # img_arr[782] represents the 782th row of the array
    # when img_arr[782] = [653, 897], that means that pixels [782,653] and [782,897] are above the limit
    # if img_arr[456] = [], that means that in the 456th row, there are no pixels that are above the limit
    img_arr = [[]] #ignore the 0th row as not having anything
    for i_row in range(1, len(img)): #range starts at 1 bc of the assumption that img[0] has no peaks
        # apply max_peaks to each row (similar to Lab2)
        img_arr.append(max_peaks(img[i_row], width=10, lower_limit=lim))

    def cluster(img_arr):
        """ The return value cluster_arr is an array where it takes in img_arr and
        sees something like this (where the number correlates to the index of img_arr 
        and the array at that index)
        287 []; 288 []; 289 []; 290 [654, 786]; 291 [656]; 292 []; 293 []; 294 []
        and groups it like so [ [290, 654], [290,786], [291,656] ] -- this is one cluster
        cluster_arr is an array of multiple 'cluster's like above
        """
        cluster_arr, k = [], 0 # index k represents the current index in img_arr
        while k < len(img_arr)-1:
            curr = img_arr[k]
            if len(curr) >= 1:
                temp_arr = [[k,elem] for elem in curr]
                j=k+1
                while j < len(img_arr) and len(img_arr[j]) >= 1:
                    for elem in img_arr[j]:
                        temp_arr.append([j,elem])
                    j+=1
                k=j-1
                cluster_arr.append(temp_arr)
            k+=1
        return cluster_arr
        
    #for index, elem in enumerate(img_arr):
    #    print index, elem
        

    clusters = cluster(img_arr)

    #for cluster in clusters:
    #    print cluster
    
    def go_thru_clusters(clusters):
        new_cluster_arr =[]
        for cluster in clusters:
            #print '*** orig cluster', cluster
            def group(cluster):
                temp_arr = [cluster[0]]
                cluster.remove(cluster[0])
                curr_i = 0
                while cluster:
                    if np.abs(cluster[curr_i][1] - temp_arr[-1][1])<20 and np.abs(cluster[curr_i][0] - temp_arr[-1][0]) <=20:
                        # first cond: deals with the diff in cols
                        # second cond: deals with the diff in rows
                        temp_arr.append(cluster[curr_i])
                        cluster.remove(cluster[curr_i])
                        #print "TEMP", temp_arr
                        #print "CLUST", cluster
                        curr_i=0
                    else:
                        curr_i+=1
                        #print curr_i
                    if len(cluster) <= curr_i:
                        break
                #print '**temp ', temp_arr
                #print '***cluster', cluster
                new_cluster_arr.append(temp_arr)
            while cluster:
                group(cluster)
        return new_cluster_arr  
    
    # new cluster arr just groups the clusters way better and applies certain restrictions
    new_cluster_arr = go_thru_clusters(clusters)
    
    # === get rid of the clusters that are length one because those are just systematic errors
    # that were picked up/or pixels that were above the limit but have no correlating pixels
    # near them
    delete_arr=[]
    for index, elem in enumerate(new_cluster_arr):
        if len(elem) == 1:
            delete_arr.append(index)
    new_cluster_arr = np.delete(new_cluster_arr, delete_arr)

    peaks = [] # each element of peaks has the syntax of [[x,y],size]
    for cluster in new_cluster_arr:
        intensity_arr = [img[coords[0]][coords[1]] for coords in cluster]
        max_index = np.argmax(intensity_arr)
        peaks.append([cluster[max_index], len(cluster)])
    return peaks

def centroid(img):
    """ 
    Purpose:
        To find the centroids of each peak, given the coordinates of each peak and 
        their corresponding size. To find the centroid, we simply compute the center
        of mass in both the x and y direction. To find x_cm and y_cm (center of mass),
        we find the x_cm of each row (or the y_cm of each col) and take the average.
        For example, if we have 3 rows and 3 columns, then we take the x_cm of the first,
        second, and third row individually and then take the average of the x_cm's (likewise 
        for y_cm)
    Input:
        fil: the uncorrected, raw file, that will be corrected
        bias: the bias 
    Output: 
        centroids_x: correlates w cols, all the x_cm's
        centroids_y: correlates w rows, all the y_cm's
        [centroids_x[i], centroids_y[i]] correlate to a centroid pixel 
    peaks has the form of [[x,y],size] 
    deciding to take into account +- size for centroiding
    """
    centroids =[]
    peaks = find_peaks(img) # find the peaks, elem of peaks is [[x,y],size]
    for pair, size in peaks:
        x_lower = 0 if pair[0]-size <= 0 else pair[0]-size
        x_upper = 1023 if pair[0]+size >=1024 else pair[0]+size
        y_lower = 0 if pair[1]-size <= 0 else pair[1]-size
        y_upper = 1023 if pair[1]+size >=1024 else  pair[1]+size
        x_vals = range(x_lower, x_upper) # just the ranging x values
        y_vals = range(y_lower, y_upper) # just the ranging y values
        x_arr =[]
        for y in y_vals:
            x_temp = []
            for x in x_vals:
                x_temp.append([x,y]) # keeps y constant for varying x; will be used for the x_cm
            x_arr.append(x_temp)
        
        y_arr=[]   
        for x in x_vals:
            y_temp = []
            for y in y_vals:
                y_temp.append([x,y]) # keeps y constant for varying y; will be used for the y_cm
            y_arr.append(y_temp)
        
        def find_cm(arr, vals):
            cm_arr=[]
            for elem in arr:
                range_img = [img[i[0],i[1]] for i in elem] # vals of img at constant y (or x), varying x (or y)
                centroid =  np.sum([a*b for a,b in zip(vals,range_img)])/np.sum(range_img) #<x> or <y>
                cm_arr.append(centroid)
            return np.sum(cm_arr)/len(cm_arr) # return the avg vals of the center of masses

        x_cm = find_cm(x_arr, x_vals)
        y_cm = find_cm(y_arr, y_vals)
        centroids.append([x_cm,y_cm])
    
    centroids = np.transpose(centroids)
    centroids_x, centroids_y = centroids[0], centroids[1]
    return centroids_x, centroids_y  # centroids_x ~ cols, centroids_y ~ rows, so when plotting in plt, 
    # need to plot plt.plot(centroids_y, centroids_x) - i think i just switched them along the way
