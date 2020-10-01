# General Stuff

## ssh-ing into baross

```bash
ssh [username ]@baross.its.carleton.edu
```


## Using screen/ tmux
(We have been using screen in class.)

See: [Cheat Sheet](http://www.dayid.org/comp/tm.html)

`^` = control key

```
| Action                     | tmux         | screen                  |
|----------------------------|--------------|-------------------------|
| start new session          | `tmux`       | `screen -S screen_name` |
| detach from current session| `^b d`       | `^a ^d`                 |
| kill current session       |              | `^a ^k`                 |
| re-attach detached session | `tmux attach`| `screen -r screen_name` |
| list sessions              | `^b s`       | `screen -ls`            |
```


## File Transfer

### Filezilla

- Open Filezilla
- Host: sftp://baross.its.carleton.edu
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

### Open via FTP with BBEdit
- In BBEdit, go to File --> Open from FTP/SFTP Server...
- Server: baross.its.carleton.edu
- Check the SFTP box
- Port: 22
- User: Carleton username
- Password: Carleton password
- Path: /Accounts/Carleton username
- Click Connect
