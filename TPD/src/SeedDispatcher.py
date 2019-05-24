import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from collections import deque 

HOST_NAME = '0.0.0.0'
PORT_NUMBER = 8080


class MyHandler(BaseHTTPRequestHandler):
    aSeed=0
    rates=[[obj]*100 for obj in range(5,100,5)]
    rates=[y for x in rates for y in x]
    indx=0
    
    def do_HEAD(self):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def do_GET(self):
        paths = {
            '/seed/new': {'status': 200},
	        '/task/new': {'status': 200},
            '/task/update': {'status': 208},
            '/baz': {'status': 404},
            '/qux': {'status': 500}
        }

        if self.path in paths:
            self.respond(paths[self.path])
        #else:
        #    self.respond({'status': 500})

    def handle_http(self, status_code, path):
        self.send_response(status_code)
        self.send_header('Content-type', 'text/html')
        self.end_headers()
        if path=='/task/new':
             if MyHandler.indx < len(MyHandler.rates):
                task=MyHandler.rates[MyHandler.indx]
                MyHandler.indx=MyHandler.indx+1
                return bytes(str(task), 'UTF-8')
             else:
                return bytes("NONE",'UTF-8')
        if path=='/seed/new':
            seed=MyHandler.aSeed+1
            MyHandler.aSeed=seed
            return bytes(str(seed),'UTF-8')

    def respond(self, opts):
        response = self.handle_http(opts['status'], self.path)
        self.wfile.write(response)

if __name__ == '__main__':
    server_class = HTTPServer
    httpd = server_class((HOST_NAME, PORT_NUMBER), MyHandler)
    print(time.asctime(), 'Server Starts - %s:%s' % (HOST_NAME, PORT_NUMBER))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        httpd.server_close()
    httpd.server_close()
    print(time.asctime(), 'Server Stops - %s:%s' % (HOST_NAME, PORT_NUMBER))
