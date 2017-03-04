#./usr/localspark-2.1.0-bin-hadoop2.7/python/
import sys
sys.path.append("/usr/local/spark/python/")
from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext, SparkSession
from pyspark.sql import *
from pyspark.sql.types import *

sc = SparkContext()
sqlCtx = SQLContext(sc)

calendar_rdd = sc.textFile("airbnb_datasets/calendar_us.csv").map(lambda x: x.split('\t')).filter(lambda line: line[0] != "listing_id").map(lambda s: (int(s[0]), s[1], s[2]))
calendar_fields=[
        StructField('listing_id', IntegerType(), False),                        #name, type, nullable
        StructField('date', StringType(), False),
        StructField('available', StringType(), False)
]

calendar_schema = StructType(calendar_fields)
calendar = sqlCtx.createDataFrame(calendar_rdd, calendar_schema)                #Create Data Frame

reviews_rdd = sc.textFile("airbnb_datasets/reviews_us.csv").map(lambda x: x.split('\t')).filter(lambda line: line[0] != "listing_id").map(lambda s: (int(s[0]), int(s[1]), s[2], int(s[3]), s[4], s[5]))
reviews_fields=[
        StructField('listing_id', IntegerType(), False),
        StructField('id', IntegerType(), False),
        StructField('date', StringType(), False),
        StructField('reviewer_id', IntegerType(), False),
        StructField('reviewer_name', StringType(), False),
        StructField('comments', StringType(), False)
]

reviews_schema = StructType(reviews_fields)
reviews = sqlCtx.createDataFrame(reviews_rdd, reviews_schema)                   #Create Data Frame


listings_rdd = sc.textFile("airbnb_datasets/listings_us.csv").map(lambda x: x.split('\t'))
listings_fields=[]

for name in range(0, len(listings_rdd.take(1)[0])):                             #Create columns to Data Frame
    temp = listings_rdd.take(1)[0][name]
    listings_fields.append(StructField(temp, StringType(), True))

listings_schema = StructType(listings_fields)
listings = sqlCtx.createDataFrame(listings_rdd, listings_schema)                #Create Data Frame

sqlCtx.registerDataFrameAsTable(calendar, "calendar")                           #Register Data Frame
sqlCtx.registerDataFrameAsTable(reviews, "reviews")
sqlCtx.registerDataFrameAsTable(listings, "listings")

print sqlCtx.sql("SELECT COUNT(DISTINCT city) AS number FROM listings").show()
sqlCtx.sql("SELECT DISTINCT city AS city FROM listings").show(200).saveAsTextFile()
