from http.server import HTTPServer, BaseHTTPRequestHandler
from collections import deque 

taskFile='E:/projects/TPD/bp2_100.txt';
tasks=deque([line.rstrip('\n') for line in open(taskFile)])

class MyHTTPRequestHandler(BaseHTTPRequestHandler):
	def do_GET(self):
		self.send_response(200)
		self.send_header('Content-type', 'text/html')
		self.end_headers()
		task=tasks.popleft()
		print(task)
		self.wfile.write(bytes(task,'utf-8'))
	
	def do_POST(self):
		self.do_GET(self)
		
	
httpd = HTTPServer(('localhost',80), MyHTTPRequestHandler)
httpd.serve_forever()