import sys

summary = "/pdb/entry/summary"
experiment = "/pdb/entry/experiment"
sifts = "/mappings"

PY3 = sys.version > '3'

if PY3:
    import urllib.request as urllib2
else:
    import urllib2

SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

def make_request(url, data):   
    request = urllib2.Request(url)

    try:
        url_file = urllib2.urlopen(request, data)
    except urllib2.HTTPError as e:
        if e.code == 404:
            print("[NOTFOUND %d] %s" % (e.code, url))
        else:
            print("[ERROR %d] %s" % (e.code, url))

        return None

    return url_file.read().decode()

def get_request(url, arg, pretty=False):
    full_url = "%s/%s/%s?pretty=%s" % (SERVER_URL, url, arg, str(pretty).lower())
    
    return make_request(full_url, None)

def post_request(url, data, pretty=False):
    full_url = "%s/%s/?pretty=%s" % (SERVER_URL, url, str(pretty).lower())
    
    if isinstance(data, (list, tuple)):
        data = ",".join(data)
    
    return make_request(full_url, data.encode())
