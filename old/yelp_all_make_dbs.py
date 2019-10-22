"""
This file is intended to create MongoDB databases for each of the yelp provided json files

author: daniel metz
"""
from pymongo import MongoClient, IndexModel, ASCENDING, DESCENDING
import ujson

client = MongoClient()
db = client.yelpdb

with open('../data/yelp_academic_dataset_business.json') as f:
    db.bus.insert_many([ujson.loads(line) for line in f], ordered=False)

with open('../data/yelp_academic_dataset_user.json') as f:
    db.user.insert_many([ujson.loads(line) for line in f], ordered=False)

with open('../data/yelp_academic_dataset_review.json') as f:
    db.review.insert_many([ujson.loads(line) for line in f], ordered=False)

# add indexes for user_id and business_id for faster access
index_uid = IndexModel([('user_id', ASCENDING)])
index_bid = IndexModel([('business_id', ASCENDING)])
db.user.create_indexes([index_uid])
db.bus.create_indexes([index_bid])
db.review.create_indexes([index_uid, index_bid])
