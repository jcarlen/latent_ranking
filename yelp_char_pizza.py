
# coding: utf-8

# In[20]:


# In[21]:

sys.path.append("/Library/Python/2.7/site-packages")
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


# In[ ]:

bus.distinct("categories")


# In[211]:

#Find all pizza places

pizza = bus.find({"categories":"Pizza"})
pizza_id = pizza.distinct('business_id')

#check
print(pizza[2])
print(pizza_id[1:3])


# In[210]:

#pizza ratings matrix
pizza_rating = np.zeros((len(pizza_id), len(pizza_id)))


# In[246]:

#find all reviews of the pizza places
rev.pizza = rev.find({'business_id':{"$in":pizza_id}}) #4808

#find all users who reviewed one of these places #2861
rev.pizza.uid = rev.pizza.distinct('user_id')

for uid in rev.pizza.uid:
    uid_reviews = list(rev.find({'user_id':uid, 'business_id':{"$in":pizza_id}}))
    x = len(uid_reviews)
    if  x > 1:
        for i in range(0,x):
            for j in range(i,x):
                diff = uid_reviews[i].get('stars') - uid_reviews[j].get('stars')
                if diff != 0:
                    I = pizza_id.index(uid_reviews[i].get('business_id'))
                    J = pizza_id.index(uid_reviews[j].get('business_id'))
                if diff > 0:   
                    pizza_rating[I,J]+=diff
                if diff < 0:
                    pizza_rating[J,I]-=diff


# In[256]:

np.savetxt("pizza_rating.csv", pizza_rating, delimiter=",")


# In[284]:

#get average rating for each pizza place
pizza_avg = range(len(pizza_id))
for k in range(0, len(pizza_id)):
    reviews = list(rev.find({'user_id':rev.pizza.uid[k]}))
    i = 0; j = 0
    for review in reviews:
        i += review['stars']
        j += 1
    pizza_avg[k] = round(float(i)/float(j),3)

print pizza_avg

np.savetxt("pizza_avg.csv", pizza_avg, delimiter=",")


# In[ ]:



