import datetime
import json
import os
import urllib

import boto3
import tweepy

TIME_UNIT = 5

s3_bucket = os.environ["DATA_S3_BUCKET"]
s3_key    = os.environ["DATA_S3_KEY"]
twitter_consumer_key    = os.environ["TWITTER_CONSUMER_KEY"]
twitter_consumer_secret = os.environ["TWITTER_CONSUMER_SECRET"]
twitter_token           = os.environ["TWITTER_TOKEN"]
twitter_token_secret    = os.environ["TWITTER_TOKEN_SECRET"]
slack_webhook_url = os.environ["SLACK_WEBHOOK_URL"]

def main(event, context):
    now = datetime.datetime.utcnow()
    for record in readDataLines(now):
        message = record[1]
        if message.startswith("#"):
            continue
        tweet(message)
        postSlack(message)
        print(f"tweet: {message}")

def readDataLines(now):
    lines1 = readS3Object().split(sep='\n')
    lines2 = [parseDataLine(line) for line in lines1]
    lines3 = [line for line in lines2 if line != None]
    lines4 = [line for line in lines3 if isTimeMatch(line, now)]
    return lines4

def isTimeMatch(record, now):
    dt = record[0]
    m0 = int(now.minute / TIME_UNIT) * TIME_UNIT
    if dt.year == now.year and dt.month == now.month and dt.day == now.day and dt.hour == now.hour:
        m1 = int(dt.minute / TIME_UNIT) * TIME_UNIT
        if m1 == m0:
            return True
    return False

def parseDataLine(line):
    line = line.strip()
    cols = line.split(sep=' ', maxsplit=1)
    try:
        dt = datetime.datetime.fromisoformat(cols[0] + "+09:00").astimezone(datetime.timezone.utc)
        msg = cols[1].strip()
    except:
        return None
    return [dt, msg]

def readS3Object():
    session = boto3.session.Session()
    s3_client = session.client("s3")
    res = s3_client.get_object(
        Bucket = s3_bucket,
        Key = s3_key,
    )
    return res['Body'].read().decode("utf-8")

def tweet(message):
    client = tweepy.Client(
        consumer_key        = twitter_consumer_key,
        consumer_secret     = twitter_consumer_secret,
        access_token        = twitter_token,
        access_token_secret = twitter_token_secret,
       )
    client.create_tweet(text = message)

def postSlack(message):
    post_data = { "text": message }
    headers = { "Content-Type": "application/json" }
    req = urllib.request.Request(slack_webhook_url, json.dumps(post_data).encode(), headers)
    with urllib.request.urlopen(req) as res:
        webhook_body = res.read()
    #print(webhook_body)

