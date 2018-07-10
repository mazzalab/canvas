var sys = require('util');
var fs = require('fs');
var path = require('path');

function walk (dir, options, callback) {
  if (!callback) {callback = options; options = {}}
  if (!callback.files) callback.files = {};
  if (!callback.pending) callback.pending = 0;
  callback.pending += 1;
  fs.stat(dir, function (err, stat) {
    if (err) return callback(err);
    callback.files[dir] = stat;
    fs.readdir(dir, function (err, files) {
      if (err) {
        if(err.code === 'EACCES' && options.ignoreUnreadableDir) return callback();
        return callback(err);
      }
      callback.pending -= 1;
      files.forEach(function (f, index) {
        f = path.join(dir, f);
        callback.pending += 1;
        fs.stat(f, function (err, stat) {
          var enoent = false
            , done = false;

          if (err) {
            if (err.code !== 'ENOENT' && (err.code !== 'EPERM' && options.ignoreNotPermitted)) {
              return callback(err);
            } else {
              enoent = true;
            }
          }
          callback.pending -= 1;
          done = callback.pending === 0;
          if (!enoent) {
            if (options.ignoreDotFiles && path.basename(f)[0] === '.') return done && callback(null, callback.files);
            if (options.filter && !options.filter(f, stat)) return done && callback(null, callback.files);
            callback.files[f] = stat;
            console.log('checking path:', f);
            if (stat.isDirectory() && !(options.ignoreDirectoryPattern && options.ignoreDirectoryPattern.test(f))) {
              console.log('walking path:', f);
              walk(f, options, callback)
            } else {
              if (stat.isDirectory()) console.log('ignoring path:', f);
            }
            done = callback.pending === 0;
            if (done) callback(null, callback.files);
          }
        })
      })
      if (callback.pending === 0) callback(null, callback.files);
    })
    if (callback.pending === 0) callback(null, callback.files);
  })

}

var options = {ignoreDirectoryPattern: new RegExp(/dist/)};

walk(".", options, function() {});
