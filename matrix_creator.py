"""
Copyright (C) 2022 Konstantin Zaytsev

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import database_sampler
import requests
from bs4 import BeautifulSoup
import time


def web_create(sample_position, number, login, password):
    database_sampler.create_file(sample_position, number)
    time_limit = 1200

    while True:
        try:
            start = time.time()
            session = requests.Session()
            # all cookies received will be stored in the session object
            auth_url = 'http://victoria.biengi.ac.ru/mahds/auth'
            auth_payload = {'login': login, 'pw': password, 'type': 'auth'}
            response = session.post(auth_url, data=auth_payload)

            mahds_url = 'http://victoria.biengi.ac.ru/mahds/main'
            file = open('bin/1000.cds', 'rb')
            mahds_payload = {'emailNotification': False, 'alignmentType': 'DNANew', 'outputType': 'Fasta'}
            response = session.post(mahds_url, data=mahds_payload, files={'seqFile': file})
            result_url = 'http://' + str(response.content)[2:-1]

            result = '\n'
            while result == '\n' and time.time() - start < time_limit:
                response = requests.get(result_url)
                soup = BeautifulSoup(response.text, features='html.parser')
                result = soup.find('div', {'class': 'mainContent'}).get_text(separator='\n')
                if result != '\n':
                    result = result.lstrip()
                    with open('bin/' + str(1000 + number) + '.txt', 'w') as file:
                        file.write(result)
                time.sleep(5)

            info = result[:300]
            class_w = (int(float(info[info.find('W22=') + 4:info.find('IM=')])))
            return class_w

        except Exception as error:
            print(error)

