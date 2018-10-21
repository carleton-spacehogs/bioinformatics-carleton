# General Stuff

## ssh-ing into liverpool / baross

Liverpool:

```bash
ssh [username]@liverpool.its.carleton.edu
```

Baross

```bash
ssh [username ]@baross.its.carleton.edu
```


## Using screen/ tmux
(We have been using screen in class.)

See: [Cheat Sheet](http://www.dayid.org/comp/tm.html)

`^` = control key

```
| Action                     | tmux         | screen     |
|----------------------------|--------------|------------|
| start new session          | `tmux`       | `screen`   |
| detach from current session| `^b d`       | `^a ^d`    |
| re-attach detached session | `tmux attach`| `screen-r` |
| list sessions              | `^b s`       | `screen-r` |
```


## File Transfer

### Filezilla

- Open Filezilla
- Host: sftp://liverpool.its.carleton.edu
- Username: Carleton username
- Password: Carleton password
- Port: 22
- Click QuickConnect

### Secure Copy

From baross to local computer:

```bash
scp [username]@baross.its.carleton.edu:/Accounts/[username]/[path of your destination directory]/[some_file.txt] ~/Desktop
```

From local computer to baross:

```bash
scp ~/Desktop/[some_file.txt] [username]@baross.its.carleton.edu:/Accounts/[username]/[path of your destination directory]
```
