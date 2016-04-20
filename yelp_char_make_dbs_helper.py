"""
make the charlotte databases (user, rev, bus)

city: charlotte

author: daniel metz
"""

import ujson
import numpy as np
from sklearn.cluster import KMeans

# load businesses
def import_businesses():
    """
        assumes the "yelp_academic_dataset_business" is in the same directory.
        
        returns the properly loaded list of json objects
        """
    with open('../data/yelp_academic_dataset_business.json') as business_file:
        businesses = [ujson.loads(line) for line in business_file]
    return businesses

# assign true_city
def assign_true_cities(businesses):
    """
        creates a new "true_city" field for each object, assigning each to one of
        ten hard-coded cities based upon their latitude and longitude
        """
    
    cluster_city_dict = make_cluster_city_dict(businesses)
    cluster_assignment = assign_clusters_cities(businesses)
    
    # add each cluster-assignment as a field for each business
    for i in range(len(businesses)):
        businesses[i]['true_city'] = cluster_city_dict[cluster_assignment[i]]
    
    # check city distribution
    city_count_dict = {}

def make_cluster_city_dict(businesses):
    """
        simply creates a hard-coded dictionary based on a deterministic KMeans
        clustering of the business cities as done in assign_clusters_cities()
        
        helper function to assign_true_cities
        """
    cluster_city_dict = {}
    cluster_city_dict[0] = 'Phoenix'
    cluster_city_dict[1] = 'Edinburgh'
    cluster_city_dict[2] = 'Montreal'
    cluster_city_dict[3] = 'Madison'
    cluster_city_dict[4] = 'Charlotte'
    cluster_city_dict[5] = 'Las Vegas'
    cluster_city_dict[6] = 'Karlsruhe'
    cluster_city_dict[7] = 'Pittsburgh'
    cluster_city_dict[8] = 'Urbana-Champaign'
    cluster_city_dict[9] = 'Waterloo'
    return cluster_city_dict

def assign_clusters_cities(businesses):
    """
        runs KMeans clustering in a deterministic fashion (uses a hard-coded seed)
        in order to assign each business to one of ten clusters
        
        helper function to assign_true_cities that simply returns a vector of
        cluster assignments
        """
    # find number of coordinate pairs
    n_samples = len(businesses)
    
    # initialize array
    coordinates = np.zeros([n_samples, 2])
    
    # fill in the coordinate array
    for sample in range(n_samples):
        coordinates[sample, 0] = businesses[sample]['latitude']
        coordinates[sample, 1] = businesses[sample]['longitude']
    
    # assign each business to a cluster
    clusters = KMeans(10, random_state=47)
    cluster_assignment = clusters.fit_predict(coordinates)
    return cluster_assignment