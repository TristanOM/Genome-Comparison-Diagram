import pymongo
from pymongo import MongoClient
import datetime
client = MongoClient('localhost' , 27017)

db = client.test_database

post = {"author": "Mike","text": "My first blog post!","tags": ["mongodb", "python", "pymongo"],"date": datetime.datetime.utcnow()}
posts = db.posts
post_id = posts.insert_one(post).insert_id
post.find_one( {u'date': datetime.datetime(...), u'text': u'My first blog post!', u'_id': ObjectId('...'), u'author': u'Mike', u'tags': [u'mongodb', u'python', u'pymongo']})
