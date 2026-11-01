// The WebGPU calls of the externs of src/gpu/gpu.goose. The C of the
// module (GPU_C, src/codegen/page.goose) calls these through
// Module.skeinGpu. This file has no solver logic: gpu.goose decides the
// buffers, the kernels of kernels.wgsl and their sequence. A handle is an
// index into the tables below. The page puts this file in the worker of
// each run, and tests/gpu/gpu.test.mjs loads it in Node.
function skeinGpu(device, wgsl) {
  const module = device.createShaderModule({ code: wgsl });
  const pipes = new Map(), bufs = [], groups = [], errors = [];
  let enc = null, pass = null;
  device.addEventListener("uncapturederror", (e) => errors.push(e.error.message));

  function pipeline(name) {
    if (!pipes.has(name)) pipes.set(name, device.createComputePipeline({ layout: "auto", compute: { module, entryPoint: name } }));
    return pipes.get(name);
  }
  function submit() {
    if (!enc) return;
    if (pass) pass.end();
    device.queue.submit([enc.finish()]);
    enc = pass = null;
  }
  return {
    errors,
    // Frees the buffers and the bind groups of the last solve.
    reset() {
      submit();
      for (const b of bufs) b.destroy();
      bufs.length = groups.length = 0;
    },
    buffer(bytes) {
      const size = Math.max(16, Math.ceil(bytes / 16) * 16);
      bufs.push(device.createBuffer({ size, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST }));
      return bufs.length - 1;
    },
    // A write goes to the queue at once, so the work recorded before it
    // goes first.
    write(b, bytes) {
      submit();
      device.queue.writeBuffer(bufs[b], 0, bytes);
    },
    // A bind group of the kernel `name`: binding pairs[2 k] is buffer pairs[2 k + 1].
    bind(name, pairs) {
      const entries = [];
      for (let k = 0; k < pairs.length; k += 2) entries.push({ binding: pairs[k], resource: { buffer: bufs[pairs[k + 1]] } });
      const pipe = pipeline(name);
      groups.push({ pipe, group: device.createBindGroup({ layout: pipe.getBindGroupLayout(0), entries }) });
      return groups.length - 1;
    },
    // A dispatch of `threads` threads, in workgroups of 64, in rows of at
    // most 65535 workgroups. The dispatches of one pass see the writes of
    // the ones before them.
    dispatch(g, threads) {
      if (!enc) enc = device.createCommandEncoder();
      if (!pass) pass = enc.beginComputePass();
      const n = Math.max(1, Math.ceil(threads / 64)), x = Math.min(n, 65535);
      pass.setPipeline(groups[g].pipe);
      pass.setBindGroup(0, groups[g].group);
      pass.dispatchWorkgroups(x, Math.ceil(n / x));
    },
    submit,
    // The first `bytes` of buffer b, after the work so far, or null after
    // a WebGPU error.
    async read(b, bytes) {
      submit();
      const size = Math.ceil(bytes / 4) * 4;
      const st = device.createBuffer({ size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
      enc = device.createCommandEncoder();
      enc.copyBufferToBuffer(bufs[b], 0, st, 0, size);
      submit();
      await st.mapAsync(GPUMapMode.READ);
      const out = new Uint8Array(st.getMappedRange(0, size)).slice(0, bytes);
      st.unmap();
      st.destroy();
      return errors.length ? null : out;
    },
  };
}
