"""
make the charlotte databases (user, rev, bus)
    
city: charlotte
    
author: daniel
"""

# load businesses
from yelp_char_make_dbs_helper import *
businesses = import_businesses()

# assign true_city
assign_true_cities(businesses)

# convert to pandas df
import pandas as pd
bus_df = pd.DataFrame(businesses)

# get only Charlotte
char_bus = bus_df[bus_df['true_city'] == 'Charlotte']

# import that into a database
from pymongo import MongoClient, IndexModel, ASCENDING, DESCENDING
client = MongoClient()
db = client.char
db.bus.insert_many([row[1].to_dict() for row in char_bus.iterrows()])

# great, so now we have the relevant businesses, now for reviews and users
# lets start with reviews, since they have user_ids and business_ids
char_bus_ids = set(db.bus.distinct('business_id'))

# read in the reviews
import ujson
with open('../data/yelp_academic_dataset_review.json') as f:
    review_list = [ujson.loads(line) for line in f]
    review_df = pd.DataFrame(review_list)

# create a list of bools corresponding to whether a given business is in
# charlotte
all_bus_ids = review_df['business_id'].ravel()
tf_list = []
for b_id in all_bus_ids:
    tf_list.append(b_id in char_bus_ids)

# insert the businesses into the database for those that are in charlotte
db.rev.insert_many([row[1].to_dict() for row in review_df[tf_list].iterrows()])

# now for users
char_user_ids = set(db.rev.distinct('user_id'))
with open('../data/yelp_academic_dataset_user.json') as f:
    user_list = [ujson.loads(line) for line in f]
    user_df = pd.DataFrame(user_list)

all_user_ids = user_df['user_id'].ravel()
tf_list = []
for u_id in all_user_ids:
    tf_list.append(u_id in char_user_ids)

db.user.insert_many([row[1].to_dict() for row in user_df[tf_list].iterrows()])

# add some indexes
index_uid = IndexModel([('user_id', ASCENDING)])
index_bid = IndexModel([('business_id', ASCENDING)])
db.user.create_indexes([index_uid])
db.bus.create_indexes([index_bid])
db.rev.create_indexes([index_uid, index_bid])
