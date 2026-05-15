import re
def test_re():
    text = '<a href="2025_10/">2025_10/</a>'
    matches = re.findall(r'href="(\d{4}_\d{2})/?', text)
    print(matches)

test_re()
