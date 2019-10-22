"""
make the R matrix for charlotte
make sub-databases for users >= 3 reviews and bus >= 5 reviews

city: Charlotte

author: Daniel Metz
"""

from pymongo import MongoClient
import numpy as np
import pandas as pd
from collections import defaultdict

# connect to Mongo
client = MongoClient()
# alias the user, rev, and bus tables
user = client.char.user
rev = client.char.rev
bus = client.char.bus

# Get review counts for each user, for each business
uids = user.distinct('user_id')
bids = bus.distinct('business_id')

uid_count = defaultdict(int)
bid_count = defaultdict(int)

for r in rev.find():
    uid_count[r['user_id']] += 1
    bid_count[r['business_id']] += 1

# get a list of user_id s and business_id s where
# the user_id is only considered ``good'' if it has >= 3 reviews
# and business_id is only considered ``good'' if it has >= 5 reviews
g_uids = []
g_bids = []

for key, value in uid_count.iteritems():
    if value >= 3:
        g_uids.append(key)

for key, value in bid_count.iteritems():
    if value >= 5:
        g_bids.append(key)

# lets sort those
g_uids.sort()
g_bids.sort()

# start the rating_diff matrix
rating_diff = np.zeros((len(g_bids), len(g_bids)))

for u in user.find():
    
    for rating in user
    
        if u['user_id'] in g_uids and r['business_id'] in g_bids:
            rating[g_uids.index(r['user_id']),
                   g_bids.index(r['business_id'])] = r['stars']

np.savetxt("rating_diff.csv", rating_diff, delimiter=",")

    # make one-to-one map of index to business_id to master_category
    # from import_businesses import *
    # from assign_master_category import *
    # businesses = import_businesses()
    # assign_master_category(businesses)
    # businesses = pd.DataFrame(businesses)
    # 
    # mapping = []
    # for b in g_bids:
    #     cat = businesses[businesses['business_id'] == b]['master_category'].item()
    #     mapping.append([g_bids.index(b), b, cat])
    # 
    # mapping = pd.DataFrame(mapping)
