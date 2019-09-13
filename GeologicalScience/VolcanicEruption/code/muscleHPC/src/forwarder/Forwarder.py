# /**
# * @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>
#
# * MUSCLE-HPC communication module
# * Copyright (C) 2016  University of Geneva, Switzerland
# *
# * MUSCLE-HPC is free software: you can redistribute it and/or
# * modify it under the terms of the GNU Affero General Public License as
# * published by the Free Software Foundation, either version 3 of the
# * License, or (at your option) any later version.
# *
# * The library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU Affero General Public License for more details.
# *
# * You should have received a copy of the GNU Affero General Public License
# * along with this program.  If not, see <http://www.gnu.org/licenses/>.
# */
# -*- coding: utf-8 -*-
import socket
import asyncore
import Queue
import struct
import uuid
import threading
lock = threading.Lock()
# code inspired from  http://code.activestate.com/recipes/483732-asynchronous-port-forwarding/

Config = {"BUFFER_SIZE": 4096}


# --------------- class HeaderInfo ----------------------------
class HeaderInfo:
    def __init__(self):
        self.operation = 0
        self.op_buffer = ''
        self.kernelNameLength_buffer = ''
        self.kernelName = ''
        self.targetEndpointLength_buffer = ''
        self.targetEndpoint = ''
        self.targetFWDLength_buffer = ''
        self.targetFWD = ''
        self.idForward = ''

    def toString(self):
        description = self.kernelName + ' -> targetEndpoint: ' + self.targetEndpoint + ' -> targetFWD: ' + self.targetFWD
        if len(self.idForward) > 0:
            description += ' idForward: ' + self.idForward
        return description

    def toBytes(self, op=1, fwd_id=''):
        operation_tmp_buf = self.op_buffer
        fwd_buf = ''
        if op == 2:  # forwarding connection
            operation_tmp_buf = struct.pack("i", op)
            fwd_buf += struct.pack("i", len(fwd_id)) + fwd_id
        return operation_tmp_buf + self.kernelNameLength_buffer + self.kernelName + self.targetEndpointLength_buffer + self.targetEndpoint + self.targetFWDLength_buffer + self.targetFWD + fwd_buf

# --------------- class utils ----------------------------


class Utils:
    def __init__(self):
        pass

    def compareUrls(self, url1, url2):
        ip1, p1 = self.convert(url1)
        ip2, p2 = self.convert(url2)
        return (ip1 == ip2) and (p1 == p2)

    def convert(self, url):
        flds = url.split(';')
        port_idx = 2
        if (len(flds) == 1):
            flds = url.split(':')
            port_idx = 1
        ip = str(flds[0])
        port = int(str(flds[len(flds) - port_idx]))
        return ip, port

# --------------- class forwarder ----------------------------


class forwarder(asyncore.dispatcher):
    def __init__(self, ip, port, remoteFwd=''):
        asyncore.dispatcher.__init__(self)
        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.set_reuse_addr()
        print ip, port
        self.myFwdAddress = ip + ';' + str(port) + ';0'
        print 'local fw=', self.myFwdAddress
        self.bind((ip, port))
        self.listen(5)
        #self.readWriters = {} # <kernelkernelName, rw>
        self.colorsRelayer = {}  # <color, rw>
        self.twinForwarder = None
        self.forwrdIdentifiers = {}  # {uniqId, rwObject}
        self.remoteForwarderRW = ''
        if remoteFwd:
            self.remoteForwarderRW = remoteFwd + ';0'
            print 'remote fw=', self.remoteForwarderRW

    def settwinForwarder(self, fwd):
        self.twinForwarder = fwd

    def handle_accept(self):
        print 'waiting'
        connection, addr = self.accept()
        with lock:
            self.processNewRequest(connection)
        #th = threading.Thread(target=self.processNewRequest, args=(connection,))
        #th.start()

    def processNewRequest(self, connection):
        print 'new connection'
        isHeader = False
        hinfo = self.readConnectionHeaderInfo(connection, isHeader)
        print hinfo.toString()
        # create readerwriter object
        if hinfo.targetEndpoint and hinfo.operation:
            print self.forwrdIdentifiers
            if (hinfo.idForward in self.forwrdIdentifiers):
                print '--> deleting request'
                #o = self.forwrdIdentifiers[hinfo.idForward]
                #initial_rw = o.remoteRW
                #o.remoteRW = None
                #o.close()
                #del self.forwrdIdentifiers[hinfo.idForward]
                #rwmanager = RemoteReaderWriter(rw=initial_rw)
                #rwmanager.connectToEndPoint(hinfo.targetEndpoint)
                readwriter = RemoteReaderWriter(conn=connection)
                rwmanager = RemoteReaderWriter(rw=readwriter)
                rwmanager.connectToEndPoint(hinfo.targetEndpoint)
            else:
                readwriter = RemoteReaderWriter(conn=connection)
                rwmanager = RemoteReaderWriter(rw=readwriter)
                status=rwmanager.connectToEndPoint(hinfo.targetEndpoint)
                if not status:  # here: hinfo.idForward is not in my list
                    if hinfo.idForward:
                        fwd_id = hinfo.idForward
                    else:
                        fwd_id = str(uuid.uuid4())
                    # add resquest id in my list-
                    self.forwrdIdentifiers.update({fwd_id: rwmanager})
                    status2 = rwmanager.connectToEndPoint(hinfo.targetFWD)
                    if status2:
                        newop = 2
                        #  forward what i have received
                        data = hinfo.toBytes(newop, fwd_id)
                        rwmanager.to_remote_queue.put(data)
                        print '--> forwarding request to:', hinfo.targetFWD, fwd_id, len(data)
                    else:
                        print '[ERROR] can not connect to', hinfo.targetFWD

                else:
                        print ' I am connected to', hinfo.targetEndpoint
        else:
            print '[Error] in request', hinfo.targetEndpoint, hinfo.operation

    def forwardRequest(self, rwmanager, hinfo, fwd_id):
        status2 = rwmanager.connectToEndPoint(hinfo.targetFWD)
        if status2:
            newop = 2
            #  forward what i have received
            data = hinfo.toBytes(newop, fwd_id)
            rwmanager.to_remote_queue.put(data)
            print '--> forwarding request to : ', hinfo.targetFWD, fwd_id, len(data)
        else:
            print '[ERROR] can not connect to', hinfo.targetFWD

    def readConnectionHeaderInfo(self, conn, isHeader):
        """Fisrt connect verify whether is to subsribe on manager or to connect
        to Relayer"""
        #if isHeader:
        #    packsize = self.readInt(conn)
        hinfo = HeaderInfo()
        hinfo.operation, hinfo.op_buffer = self.readInt(conn)
        length, hinfo.kernelNameLength_buffer = self.readInt(conn)
        hinfo.kernelName = self.readBytes(conn, length)
        length, hinfo.targetEndpointLength_buffer = self.readInt(conn)
        hinfo.targetEndpoint = self.readBytes(conn, length)
        length, hinfo.targetFWDLength_buffer = self.readInt(conn)
        hinfo.targetFWD = self.readBytes(conn, length)
        if hinfo.operation == 2:  # i.e forwarding
            length, r = self.readInt(conn)
            hinfo.idForward = self.readBytes(conn, length)
        return hinfo

    def readInt(self, conn):
        read = conn.recv(4)
        value = int(struct.unpack('i', read[0:4])[0])
        return value, read

    def readBytes(self, conn, count):
        bytesleft = count
        recvBuffer = ''
        while True:
            buff = conn.recv(bytesleft)
            rc = len(buff)
            recvBuffer += buff
            while rc <= 0:
                buff = conn.recv(bytesleft)
                rc = len(buff)
                recvBuffer += buff
                if (rc == 0):
                        break
            bytesleft -= rc
            if(not bytesleft):
                break
        return recvBuffer

# --------------- class RemoteReaderWriter ----------------------------


class RemoteReaderWriter(asyncore.dispatcher):

    def __init__(self, conn=None, rw=None):
        self.isConnected = False
        if conn:
            asyncore.dispatcher.__init__(self, conn)
            self.isConnected = True
        else:
            asyncore.dispatcher.__init__(self)
        self.to_remote_buffer = ''
        self.to_remote_queue = Queue.Queue()  # FIFO
        self.remoteRW = rw
        if rw:
            rw.remoteRW = self

    def connectToEndPoint(self, targetEndpoint):
         #[0:len(targetEndpoint)-1] # remove the last char \00x
        data = targetEndpoint.rstrip('\x00')
        """ connect to the target"""
        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.settimeout(3)   # 5 seconds
        flds = data.split(';')
        port_idx = 2
        if (len(flds) == 1):
            flds = data.split(':')
            #print 'data=', data
            port_idx = 1
        #print 'flds:', flds
        remoteport = int(str(flds[len(flds) - port_idx]))
        l = len(flds) - port_idx
        for i in range(0, l):
            remoteaddr = flds[i]
            ttr = 5
            while ttr > 0:
                try:
                    self.connect((remoteaddr, remoteport))
                    self.isConnected = True
                    return True
                except:
                    print 'can not connect to: ', remoteaddr, ':', remoteport
                ttr -= 1
        return False

    def handle_connect(self):
        pass

    def handle_read(self):
        if self.remoteRW:
            read = self.recv(Config["BUFFER_SIZE"])
            self.remoteRW.to_remote_queue.put(read)

    def writable(self):
        return ((not self.to_remote_queue.empty()) or (len(self.to_remote_buffer) > 0))

    def handle_write(self):
        if self.isConnected:
            if(len(self.to_remote_buffer) == 0):
                self.to_remote_buffer = self.to_remote_queue.get()
            sent = self.send(self.to_remote_buffer)
            self.to_remote_buffer = self.to_remote_buffer[sent:]

    def handle_close(self):
        print 'RemoteReaderWriter closed'
        self.close()

# --------------- MAIN  ----------------------------
if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    parser.add_option(
        '-l', '--local-ip',
        dest='local_ip', default='0.0.0.0',
        help='Local IP address to bind to')
    parser.add_option(
        '-p', '--local-port',
        type='int', dest='local_port', default=7000,
        help='Local port to bind to')
    parser.add_option(
        '-r', '--remote-ip', dest='remote_ip', default='0.0.0.0',
        help='Local IP address to bind to')
    parser.add_option(
        '-P', '--remote-port',
        type='int', dest='remote_port', default=9000,
        help='Remote port to bind to')
    parser.add_option(
        '-f', '--remote-fw',
        dest='remote_fw', default='',
        help='Remote IP address of Forwarder to connect to')

    options, args = parser.parse_args()

    fw1 = forwarder(options.local_ip, options.local_port, options.remote_fw)
    #fw2=forwarder(options.remote_ip,options.remote_port)
    #fw1.settwinForwarder(fw2)
    #fw2.settwinForwarder(fw1)

    asyncore.loop()
